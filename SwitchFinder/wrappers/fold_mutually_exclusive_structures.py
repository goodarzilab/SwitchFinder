import argparse
import os
import subprocess
from subprocess import PIPE, run, Popen, call, list2cmdline
import numpy as np
import re
import sys
from datetime import datetime
import pickle
import multiprocessing

import sys
sys.path.append('/switchfinder/')

import SwitchFinder.folding_api as folding_api


def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="the fasta file with sequences of multiple fragments", type=str)
    parser.add_argument("--shape_dict", help="dictionary with SHAPE/DMS values", type=str)
    parser.add_argument("-s", help="file with MI scores for multiple fragments", type=str)
    parser.add_argument("-o", help="output filename", type=str)
    parser.add_argument("--num_processes", help="number of processes to run in parralel", type=int)
    parser.add_argument("--temp_files_folder", help="for files like RNA_name*MI.txt", type=str)
    parser.add_argument("--RNAstructure_path", help="path to RNAstructure program", type=str)
    parser.set_defaults(
                        shape_dict=None,
                        f='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/UTRs/test_inputs_outputs/VEGFA_3_utr_shuffled_long.fa',
                        s='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/UTRs/test_inputs_outputs/VEGFA_10_shuffling_scores_unordered.txt',
                        temp_files_folder='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/folding_landscapes/temp_files',
                        o='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/folding_landscapes/test_outputs/ENST00000372077_part_4_original_energies.txt',
                        RNAstructure_path='/avicenna/khorms/programs/RNAstructure',
                        num_processes=10
                        )
    args = parser.parse_args(raw_args)
    return args

def read_arguments(raw_args):
    args = handler(raw_args)
    fasta_filename = args.f
    scores_filename = args.s
    output_filename = args.o
    if args.shape_dict is None:
        input_shape_dict = None
    else:
        input_shape_dict = pickle.load(open(args.shape_dict, 'rb'))

    global RNAstructure_path
    global temp_files_folder
    global NUM_PROCESSES

    NUM_PROCESSES = args.num_processes
    temp_files_folder = args.temp_files_folder
    RNAstructure_path = os.path.join(args.RNAstructure_path, 'exe')
    RNAstructure_parameters_data_path = os.path.join(args.RNAstructure_path, 'data_tables')
    # RNAstructure won't work without DATAPATH environmental variable set to the right folder
    # I can't set it with subprocess.Popen for reasons described here: https://stackoverflow.com/questions/3929319/subprocess-module-errors-with-export-in-python-on-linux
    os.putenv("DATAPATH", RNAstructure_parameters_data_path)

    return fasta_filename, input_shape_dict, scores_filename, output_filename


def write_files_for_folding(fragment_name):
    loops_dict_fr = folding_api.get_two_major_conflicting_loops(init_pars_dict_loc[fragment_name])
    constraints_filenames = folding_api.write_all_constraints_files_RNAstructure(loops_dict_fr, fragment_name,
                                                                                 temp_files_folder)
    sequence_filename = folding_api.write_fragment_sequence_file(fragment_name, fragments_seq_dict_loc[fragment_name],
                                                     temp_files_folder)
    if not shape_dict is None:
        shape_filename = folding_api.write_shape_constraints_RNAstructure(fragment_name, shape_dict[fragment_name],
                                                                          temp_files_folder)
    else:
        shape_filename = None

    return constraints_filenames, sequence_filename, shape_filename

def worker_for_parallel_implementation(fragment_name):
    constraints_filenames, sequence_filename, shape_filename = write_files_for_folding(fragment_name)
    filenames_to_remove = constraints_filenames + [sequence_filename]
    dotbrackets_list = []
    ct_files = []

    for constr_filename in constraints_filenames:
        pfs_filename = constr_filename.replace(".txt", "_partition.pfs")
        ct_filename = constr_filename.replace(".txt", "_MEA.ct")
        partition_command = os.path.join(RNAstructure_path, "partition")
        MEA_command = os.path.join(RNAstructure_path, "MaxExpect")
        partition_arguments = [partition_command, sequence_filename, pfs_filename, '--constraint', constr_filename]
        MEA_arguments = [MEA_command, pfs_filename, ct_filename]

        if not shape_filename is None:
            partition_arguments.append('--SHAPE')
            partition_arguments.append(shape_filename)
            filenames_to_remove.append(shape_filename)

        partition_result = run(args=partition_arguments,  stdout=PIPE, stderr=PIPE, cwd=temp_files_folder)
        MEA_result = run(args=MEA_arguments,  stdout=PIPE, stderr=PIPE, cwd=temp_files_folder)
        MEA_bracket_notation, dot_filename = folding_api.ctfile_to_dotbracket(ct_filename,
                                                                              RNAstructure_path)
        dotbrackets_list.append(MEA_bracket_notation)
        ct_files.append(ct_filename)
        filenames_to_remove += [pfs_filename, ct_filename, dot_filename]

    common_unwound_structure_file = folding_api.write_common_unwound_structure(fragment_name,
                                                             temp_files_folder,
                                                             ct_files)
    ct_files.append(common_unwound_structure_file)
    common_bracket_notation, common_dot_filename = folding_api.ctfile_to_dotbracket(common_unwound_structure_file,
                                                                            RNAstructure_path)
    filenames_to_remove += [common_unwound_structure_file, common_dot_filename]

    consensus_structures_dict = {}
    consensus_structures_dict['major'] = dotbrackets_list[0]
    consensus_structures_dict['second'] = dotbrackets_list[1]
    consensus_structures_dict['common'] = common_bracket_notation
    free_energies_dict, energy_filenames = folding_api.calculate_free_energies(ct_files,
                                                                               RNAstructure_path)
    filenames_to_remove += energy_filenames
    string_to_write = folding_api.make_sting_to_write(fragment_name,
                                                      fragments_seq_dict_loc[fragment_name],
                                                      free_energies_dict,
                                                      consensus_structures_dict)
    for fn in sorted(list(set(filenames_to_remove))):
        os.remove(fn)
    return string_to_write




def get_neccesary_dictionaries(scores_filename, input_shape_dict, fasta_filename):
    global init_pars_dict_loc
    global fragments_seq_dict_loc
    global shape_dict

    init_pars_dict_loc = folding_api.conflicting_loops_scores_output_parser(scores_filename)
    fragments_seq_dict_loc = folding_api.get_all_the_fragments_from_fasta_file(fasta_filename)
    if input_shape_dict is None:
        shape_dict = None
    else:
        shape_dict = input_shape_dict.copy()


def main(raw_args = None):
    fasta_filename, input_shape_dict, scores_filename, output_filename = read_arguments(raw_args)
    get_neccesary_dictionaries(scores_filename, input_shape_dict, fasta_filename)
    pool = multiprocessing.Pool(NUM_PROCESSES)
    fragment_names_to_iterate_through = list(init_pars_dict_loc.keys())

    with open(output_filename, 'w') as wf:
        for string_to_write in pool.imap_unordered(worker_for_parallel_implementation, fragment_names_to_iterate_through):
            wf.write(string_to_write)
            wf.flush()


if __name__ == "__main__":
    main()