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


import SwitchFinder.perturbations as perturbations
import SwitchFinder.utils as utils
import SwitchFinder.IO as IO
import SwitchFinder.folding_api as folding_api


def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="the fasta file with sequences of multiple fragments", type=str)
    parser.add_argument("-s", help="file with MI scores for multiple fragments", type=str)
    parser.add_argument("-o", help="output folder", type=str)
    parser.add_argument("--const_seq_left", help="", type=str)
    parser.add_argument("--const_seq_right", help="", type=str)
    parser.add_argument("--min_fraction_unpaired_second_loop", help="", type=float)
    parser.add_argument("--max_stretch_base_pairs_full_seq", help="", type=int)
    parser.add_argument("--n_iterations", help="", type=int)
    parser.add_argument("--n_loops_desired", help="", type=int)
    parser.add_argument("--num_processes", help="number of processes to run in parralel", type=int)
    parser.set_defaults(
                        f='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/UTRs/test_inputs_outputs/VEGFA_3_utr_shuffled_long.fa',
                        s='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/UTRs/test_inputs_outputs/VEGFA_10_shuffling_scores_unordered.txt',
                        o='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/folding_landscapes/test_outputs/ENST00000372077_part_4_original_energies.txt',
                        const_seq_left="AAAAAAAAAAAAAAA",
                        const_seq_right="AAAAAAAAAAAAAAA",
                        min_fraction_unpaired_second_loop = 0.6,
                        max_stretch_base_pairs_full_seq = 3,
                        n_iterations = 10000,
                        n_loops_desired = 10,
                        num_processes = 10,
                        how_often_to_print = 10,
                        )
    args = parser.parse_args(raw_args)
    return args

def read_arguments(raw_args = None):
    args = handler(raw_args)
    fasta_filename = args.f
    scores_filename = args.s
    output_folder = args.o
    n_processes = args.num_processes
    how_often_to_print = args.how_often_to_print

    short_filename = fasta_filename.split('/')[-1].split('.')[0]
    output_filename = os.path.join(output_folder, 'perturbations.txt')
    output_timing_filename = os.path.join(output_folder, short_filename + 'timing.txt')

    return fasta_filename, scores_filename, \
           output_filename, output_timing_filename, \
           n_processes, how_often_to_print, args


def worker_for_parallel_implementation(fr_name):
    curr_conflict_dict = fragments_conflicts_dict[fr_name]
    curr_seq_string = fragments_seq_dict[fr_name]
    all_perturbations_dict = perturbations.all_perturbations_wrapper_one_fragment(
                            curr_conflict_dict,
                            curr_seq_string,
                            const_seq_left,
                            const_seq_right,
                            min_fraction_unpaired_second_loop = min_fraction_unpaired_second_loop,
                            max_stretch_base_pairs_full_seq = max_stretch_base_pairs_full_seq,
                            n_iterations = n_iterations,
                            n_loops_desired = n_loops_desired
                        )
    string_to_write = IO.write_possible_perturbations_dict(fr_name, all_perturbations_dict)
    return string_to_write


def set_global_variables(scores_filename, fasta_filename, args):
    # a potential solution https://stackoverflow.com/questions/10152377/alternative-use-patterns-for-python-multiprocessing-avoiding-proliferation-of-gl
    # implement in the next script
    global fragments_conflicts_dict
    global fragments_seq_dict

    init_pars_dict_loc = folding_api.conflicting_loops_scores_output_parser(scores_filename)
    fragments_conflicts_dict = utils.get_confl_loops_dict(init_pars_dict_loc)
    fragments_seq_dict = folding_api.get_all_the_fragments_from_fasta_file(fasta_filename)

    global const_seq_left
    global const_seq_right
    global min_fraction_unpaired_second_loop
    global max_stretch_base_pairs_full_seq
    global n_iterations
    global n_loops_desired

    const_seq_left = args.const_seq_left.upper()
    const_seq_right = args.const_seq_right.upper()
    min_fraction_unpaired_second_loop = args.min_fraction_unpaired_second_loop
    max_stretch_base_pairs_full_seq = args.max_stretch_base_pairs_full_seq
    n_iterations = args.n_iterations
    n_loops_desired = args.n_loops_desired


def mp_handler(fasta_filename, scores_filename,
               output_filename, output_timing_filename,
               n_processes, how_often_to_print, args):
    set_global_variables(scores_filename, fasta_filename, args)
    pool = multiprocessing.Pool(n_processes)
    fragment_names_to_iterate_through = list(fragments_conflicts_dict.keys())
    counter = 0

    with open(output_filename, 'w') as wf:
        with open(output_timing_filename, 'w') as write_timing:
            print("Started")
            write_timing.write("Beginning the calculation:\t" +
                               datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
            for string_to_write in pool.imap_unordered(worker_for_parallel_implementation, fragment_names_to_iterate_through):
                print('Done: fragment number ', counter)
                counter += 1
                wf.write(string_to_write)
                wf.flush()
                if counter % how_often_to_print == 0:
                    write_timing.write(("%d fragments have been processed:\t" % (counter)) +
                                       datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
                    write_timing.flush()



def main(raw_args = None):
    fasta_filename, scores_filename, \
    output_filename, output_timing_filename,\
    n_processes, how_often_to_print, args = read_arguments(raw_args)

    mp_handler(fasta_filename, scores_filename,
               output_filename, output_timing_filename,
               n_processes, how_often_to_print, args)



if __name__ == "__main__":
    main()
