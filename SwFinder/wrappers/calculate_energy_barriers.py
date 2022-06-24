import argparse
import os
import subprocess
import shutil
import numpy as np
import re
import math
import pandas as pd
from datetime import datetime
import multiprocessing

CONSTANT_FOR_FAILING = 1000000

def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dotbracket", help="the output file for ", type=str)
    parser.add_argument("-o", help="output filename", type=str)
    parser.add_argument("--num_processes", help="number of processes to run in parralel", type=int)
    parser.add_argument("--temp_files_folder", help="for files like RNA_name*MI.txt", type=str)
    parser.add_argument("--path_rnapathfinder", help="", type=str)
    parser.add_argument("--how_often_print", help="print when each N-th fragment is done", type=int)
    parser.set_defaults(
                        dotbracket='/Users/student/Documents/hani/SNIP_switches/MIBP_python/UTRs/test_inputs_outputs/example_energy_file.txt',
                        o='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/folding_landscapes/test_outputs/ENST00000372077_part_4_original_energies.txt',
                        path_rnapathfinder='/avicenna/khorms/programs/RNApathfinder/srcTABU/get_barrier',
                        num_processes=10,
                        how_often_print=100
                        )
    args = parser.parse_args()
    return args


def read_arguments():
    args = handler()
    dotbracket_filename = args.dotbracket
    output_filename = args.o

    global NUM_PROCESSES
    global temp_files_folder
    global PATH_THAPATHFINDER
    global HOW_OFTEN_TO_PRINT

    temp_files_folder = args.temp_files_folder
    NUM_PROCESSES = args.num_processes
    PATH_THAPATHFINDER = args.path_rnapathfinder
    HOW_OFTEN_TO_PRINT = args.how_often_print

    return dotbracket_filename, output_filename


def read_dotbracket_file(dotbracket_filename):
    fragments_consensus_structures_dict = {}
    fragments_sequence_dict = {}
    with open(dotbracket_filename, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n$$$\n')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            splitted_entry = entry.split('\n')
            fragment_name = splitted_entry[0]
            fragment_sequence = splitted_entry[1]

            fragments_consensus_structures_dict[fragment_name] = {}
            fragments_sequence_dict[fragment_name] = fragment_sequence

            fragments_consensus_structures_dict[fragment_name]['major_loop'] = splitted_entry[3]
            fragments_consensus_structures_dict[fragment_name]['second_loop'] = splitted_entry[5]
            fragments_consensus_structures_dict[fragment_name]['common'] = splitted_entry[7]

    return fragments_consensus_structures_dict, fragments_sequence_dict


def get_neccesary_dictionaries(dotbracket_filename):
    global fragments_consensus_structures_dict
    global fragments_sequence_dict

    fragments_consensus_structures_dict, fragments_sequence_dict = read_dotbracket_file(dotbracket_filename)


def compile_global_patterns():
    global energy_pattern_RNAeval
    energy_pattern_RNAeval = re.compile("\(\s*\-*\d+\.\d+\)")


def write_fragment_sequence_file(fragment_name, current_structure, sequence, current_filename):
    with open(current_filename, 'w') as wf:
        wf.write(">%s\n%s\n%s\n" % (fragment_name, sequence, current_structure))


def parse_RNAeval_output(curr_filename):
    with open(curr_filename, 'r') as rf:
        big_line = rf.read().split('\n')
        line_of_interest = big_line[2]
        res = energy_pattern_RNAeval.search(line_of_interest)
        if not res:
            # print(line_of_interest)
            print("Error")
            return
        else:
            energy = float(res.group()[1:-1])

    return energy


def calculate_structure_free_energy_RNAeval(fragment_name, current_structure, sequence, loop_type):
    seq_str_filename = os.path.join(temp_files_folder, "%s_%s_sequence.fa" % (fragment_name, loop_type))
    write_fragment_sequence_file(fragment_name, current_structure, sequence, seq_str_filename)
    fold_out_file = seq_str_filename.replace("_sequence.fa", "_RNAeval_output.txt")
    command_to_run = "RNAeval < %s" % (seq_str_filename)
    with open(os.path.join(fold_out_file), 'w') as fout:
        subprocess.call(command_to_run, stdout=fout, shell=True)
    total_energy = parse_RNAeval_output(fold_out_file)
    os.remove(seq_str_filename)
    os.remove(fold_out_file)
    return total_energy


def parse_RNApathfinder_output(rnapath_out_filename):
    with open(rnapath_out_filename, 'r') as rf:
        for line in rf:
            if line.startswith("barrier is "):
                splitted_line = line.split(' -- ')
                barrier_str = splitted_line[0].replace("barrier is ", "")
                barrier_float = float(barrier_str)
                return barrier_float
    return CONSTANT_FOR_FAILING


def run_RNApathfinder(current_consensus_dict, free_energies_dict, fragment_name, fragment_sequence):
    free_energies_dict_barriers = free_energies_dict.copy()
    energy_boundary = max((free_energies_dict['common'] - free_energies_dict['major_loop']),
                        (free_energies_dict['common'] - free_energies_dict['second_loop']))
    energy_boundary += 1

    RNApathfinder_command = "%s %s '%s' '%s' 1 10 70 %d" % (PATH_THAPATHFINDER, fragment_sequence, 
                                                            current_consensus_dict['major_loop'],
                                                            current_consensus_dict['second_loop'],
                                                            energy_boundary
                                                            )
    rnapath_out_filename = os.path.join(temp_files_folder, "%s_RNApathfinder_output.txt" % (fragment_name))
    with open(rnapath_out_filename, 'w') as fout:
        subprocess.call(RNApathfinder_command, stdout=fout, shell=True)

    barrier_energy = parse_RNApathfinder_output(rnapath_out_filename)
    free_energies_dict_barriers['common'] = barrier_energy
    os.remove(rnapath_out_filename)
    return free_energies_dict_barriers


def calculate_all_energies_RNAeval(fragment_name, current_consensus_dict, fragment_sequence):
    free_energies_dict = {}
    major_loop_free_energy = calculate_structure_free_energy_RNAeval(fragment_name, 
                                                                        current_consensus_dict['major_loop'], 
                                                                        fragment_sequence, 'major_loop')
    second_loop_free_energy = calculate_structure_free_energy_RNAeval(fragment_name, 
                                                                        current_consensus_dict['second_loop'], 
                                                                        fragment_sequence, 'second_loop')
    common_free_energy = calculate_structure_free_energy_RNAeval(fragment_name, 
                                                                        current_consensus_dict['common'], 
                                                                        fragment_sequence, 'common')
    free_energies_dict['major_loop'] = major_loop_free_energy
    free_energies_dict['second_loop'] = second_loop_free_energy
    free_energies_dict['common'] = common_free_energy
    return free_energies_dict


def make_sting_to_write(fragment_name, fragment_sequence, free_energies_dict,
                        consensus_structures_dict):

    string_to_write = '%s\n%s\n' % (fragment_name, fragment_sequence)
    string_to_write += 'major_loop\n%s\nsecond_loop\n%s\ncommon\n%s\n' % (consensus_structures_dict['major_loop'], consensus_structures_dict['second_loop'], consensus_structures_dict['common'])
    string_to_write += 'major_loop: %f, second_loop: %f, common: %f\n' % (free_energies_dict['major_loop'], free_energies_dict['second_loop'], free_energies_dict['common'])
    string_to_write += ''
    string_to_write += '$$$\n'
    return string_to_write


def worker_for_parallel_implementation(fragment_name):
    fragment_sequence = fragments_sequence_dict[fragment_name]

    
    free_energies_dict = calculate_all_energies_RNAeval(fragment_name,
                                                        fragments_consensus_structures_dict[fragment_name], 
                                                        fragment_sequence)

    free_energies_dict_barriers = run_RNApathfinder(fragments_consensus_structures_dict[fragment_name], 
                                                                    free_energies_dict, fragment_name, 
                                                                    fragment_sequence)
    
    string_to_write = make_sting_to_write(fragment_name, fragment_sequence, free_energies_dict_barriers,
                                          fragments_consensus_structures_dict[fragment_name])
    return string_to_write


def main():
    dotbracket_filename, output_filename = read_arguments()
    get_neccesary_dictionaries(dotbracket_filename)
    compile_global_patterns()
    pool = multiprocessing.Pool(NUM_PROCESSES)
    fragment_names_to_iterate_through = list(fragments_consensus_structures_dict.keys())

    counter = 0
    with open(output_filename, 'w') as wf:
        print("Beginning the calculation:\t" + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
        for string_to_write in pool.imap_unordered(worker_for_parallel_implementation, fragment_names_to_iterate_through):
            counter += 1
            wf.write(string_to_write)
            wf.flush()
            if counter % HOW_OFTEN_TO_PRINT == 0:
                print("%d fragments have been processed:\t" % (counter) + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')



if __name__ == "__main__":
    main()