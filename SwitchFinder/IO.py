import os
import re
import pandas as pd

import sys
sys.path.append('/switchfinder/')

from SwitchFinder.classes import fragment, conformation, fragment_collection
import SwitchFinder.utils as sw_utils
import SwitchFinder.perturbations as perturbations

def entry_to_fragment_object(entry):
    current_fragment = fragment()

    splitted_entry = entry.split('\n')
    current_fragment.name = splitted_entry[0]
    current_fragment.sequence = splitted_entry[1]
    current_fragment.major_conf.string = splitted_entry[3]
    current_fragment.second_conf.string = splitted_entry[5]
    current_fragment.common_conf.string = splitted_entry[7]

    return current_fragment


def read_dotbracket_file(dotbracket_filename):
    fragments_consensus_structures = fragment_collection()

    with open(dotbracket_filename, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n$$$\n')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            current_fragment = entry_to_fragment_object(entry)
            fragments_consensus_structures.body_dict[current_fragment.name] = current_fragment

    return fragments_consensus_structures


def read_fasta(infile):
    tr_dict_loc = {}
    with open(infile, 'r') as f:
        split_string = f.read().split('>')
        for entry in split_string:
            if entry == '':
                continue
            seq_start = entry.find('\n')
            annotation = entry[:seq_start]
            sequence = entry[seq_start + 1:].replace('\n', '')
            tr_dict_loc[annotation] = sequence

    return tr_dict_loc


def read_folded_fasta(infile):
    folded_dict_loc = {}
    with open(infile, 'r') as f:
        split_string = f.read().split('>')
        for entry in split_string:
            if entry == '':
                continue
            split_entry = entry.split('\n')
            assert len(split_entry) == 4, "Wrong formatting! There should be 3 lines per each fragment (and an empty line)"
            annotation = split_entry[0]
            sequence = split_entry[1]
            folding = split_entry[2]
            folded_dict_loc[annotation] = {
                "sequence" : sequence,
                "dotbracket" : folding
            }
    return folded_dict_loc


def write_fasta(inp_dict, filename):
    with open(filename, 'w') as wf:
        for i in sorted(list(inp_dict.keys())):
            string_to_write = ">%s\n%s\n" % (i, inp_dict[i])
            wf.write(string_to_write)


def write_possible_perturbations_dict(fr_name, perturb_dict,
                                      inserts_as_strings = False):
    seq_np = perturb_dict["sequence_np"]
    seq_str = sw_utils.array_to_string(seq_np)
    conflict_dict = perturb_dict["conflict"]

    major_loop, second_loop = perturbations.initialize_two_loops(
                                    conflict_dict['major'],
                                    conflict_dict['second'],
                                    seq_np)

    major_string = major_loop.get_dot_bracket_string(seq_np)
    second_string = second_loop.get_dot_bracket_string(seq_np)

    string_to_write = ""
    string_to_write += "%s\n" % (fr_name)
    string_to_write += "%s\n" % (seq_str)
    string_to_write += "%s | %s\n" % (str(conflict_dict['major']["forward"]), str(conflict_dict['major']["backward"]))
    string_to_write += "%s | %s\n" % (str(conflict_dict['second']["forward"]), str(conflict_dict['second']["backward"]))
    string_to_write += "%s\n" % (major_string)
    string_to_write += "%s\n" % (second_string)
    if not inserts_as_strings:
        string_to_write += "major_strengthen\n"
        string_to_write += "%s\n" % "\n".join(perturbations.list_of_loops_to_sequences(perturb_dict["major_strengthen"], seq_np))
        string_to_write += "second_strengthen\n"
        string_to_write += "%s\n" % "\n".join(perturbations.list_of_loops_to_sequences(perturb_dict["second_strengthen"], seq_np))
        string_to_write += "second_weaken\n"
        string_to_write += "%s\n" % "\n".join(perturbations.list_of_loops_to_sequences(perturb_dict["second_weaken"], seq_np))
        string_to_write += "major_weaken\n"
        string_to_write += "%s\n" % "\n".join(perturbations.list_of_loops_to_sequences(perturb_dict["major_weaken"], seq_np))
    else:
        string_to_write += "major_strengthen\n"
        string_to_write += "%s\n" % "\n".join(perturb_dict["major_strengthen"])
        string_to_write += "second_strengthen\n"
        string_to_write += "%s\n" % "\n".join(perturb_dict["second_strengthen"])
        string_to_write += "second_weaken\n"
        string_to_write += "%s\n" % "\n".join(perturb_dict["second_weaken"])
        string_to_write += "major_weaken\n"
        string_to_write += "%s\n" % "\n".join(perturb_dict["major_weaken"])
    string_to_write += "$$$\n"

    return string_to_write


def parse_one_loop_line(line):
    numbers_in_string = re.findall(r'\d+\.*\d*', line)
    numbers_in_string = [float(x) for x in numbers_in_string]
    curr_dict = {
        "forward": (numbers_in_string[0], numbers_in_string[1]),
        "backward": (numbers_in_string[2], numbers_in_string[3]),
        "average_MI" : 0
    }
    return curr_dict



def read_possible_perturbation_file(inp_filename):
    perturbations_dict = {}
    with open(inp_filename, 'r') as rf:
        big_string = rf.read()
    big_array = big_string.split("$$$")
    for string_item in big_array:
        string_item_splitted = string_item.split('\n')
        lines = [x for x in string_item_splitted if x != ""]
        if len(lines) == 0:
            continue
        fr_name = lines[0]

        # get the sequence
        sequence_str = lines[1]
        sequence_np = sw_utils.string_to_array(sequence_str)

        # get the conflict dict
        major_loop_coordinates = parse_one_loop_line(lines[2])
        second_loop_coordinates = parse_one_loop_line(lines[3])
        conflict_dict = {"major" : major_loop_coordinates,
                         "second" : second_loop_coordinates}

        # get the dotbracket structures - just in case
        dotbracket_major = lines[4]
        dotbracket_second = lines[5]

        # get the perturbations dict
        major_strengthen_index = lines.index('major_strengthen')
        second_strengthen_index = lines.index('second_strengthen')
        second_weaken_index = lines.index('second_weaken')
        major_weaken_index = lines.index('major_weaken')

        major_strengthen_lines = lines[major_strengthen_index + 1 :
                                 second_strengthen_index]
        second_strengthen_lines = lines[second_strengthen_index + 1 :
                                 second_weaken_index]
        second_weaken_lines = lines[second_weaken_index + 1 :
                                 major_weaken_index]
        # in case there are no major_weaken lines and we are getting out of the list
        major_weaken_lines = lines[min(major_weaken_index + 1, len(lines)) :
                                 ]
        perturbations_dict[fr_name] = {
            "conflict" : conflict_dict,
            "sequence_np" : sequence_np,
            "dotbracket_major" : dotbracket_major,
            "dotbracket_second" : dotbracket_second,
            "major_strengthen" : major_strengthen_lines,
            "second_strengthen" : second_strengthen_lines,
            "second_weaken" : second_weaken_lines,
            "major_weaken" : major_weaken_lines,
        }

    return perturbations_dict


def single_perturbation_list_to_strings(fr_name,
                                        perturbation_type,
                                        perturbations_list):
    sequences_list = []
    for i in range(len(perturbations_list)):
        current_name = "%s|%s|%s" % (fr_name, perturbation_type, i)
        current_string = ">%s\n%s\n" % (current_name, perturbations_list[i])
        sequences_list.append(current_string)
    return sequences_list


def write_out_perturbed_sequences(perturb_dict):
    pert_types_list = ["major_strengthen",
                       "second_strengthen",
                       "second_weaken",
                       "major_weaken"]
    strings_to_write_list = []
    for fr_name in sorted(list(perturb_dict.keys())):
        for p_type in pert_types_list:
            curr_list = single_perturbation_list_to_strings(
                        fr_name,
                        p_type,
                        perturb_dict[fr_name][p_type])
            strings_to_write_list += curr_list
    string_to_write = "".join(strings_to_write_list)
    return string_to_write


def read_selection_dataframe(inp_filename):
    curr_df = pd.read_csv(inp_filename, sep='\t')
    return curr_df

def create_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)