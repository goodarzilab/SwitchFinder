import sys
import numpy as np
import copy
import random
import pandas as pd
import io
import time
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.genmod.families import Binomial
import statsmodels.formula.api as smf
import statsmodels.api as sm
import itertools
from scipy.stats import pearsonr
# from pandarallel import pandarallelmutated_nt_content_difference

import sys
sys.path.append('/switchfinder/')

import SwitchFinder.utils as sw_utils
import SwitchFinder.IO as sw_io
import SwitchFinder.glob_vars
from SwitchFinder.classes import loop, base_pair_array
import SwitchFinder.folding_api



def convert_perturbation_dict_to_df(perturbations_dict):
    potential_sequences_found_dict = {}
    for fr in perturbations_dict:
        potential_sequences_found_dict[fr] = {}
        potential_sequences_found_dict[fr]["major_strengthen"] = len(perturbations_dict[fr]["major_strengthen"]) > 0
        potential_sequences_found_dict[fr]["second_strengthen"] = len(perturbations_dict[fr]["second_strengthen"]) > 0
        potential_sequences_found_dict[fr]["second_weaken"] = len(perturbations_dict[fr]["second_weaken"]) > 0
        potential_sequences_found_dict[fr]["major_weaken"] = len(perturbations_dict[fr]["major_weaken"]) > 0

        first_loop, second_loop = initialize_two_loops(
                                            perturbations_dict[fr]['conflict']['major'],
                                            perturbations_dict[fr]['conflict']['second'],
                                            perturbations_dict[fr]["sequence_np"])
        potential_sequences_found_dict[fr]["first_loop_length"] = first_loop.length
        potential_sequences_found_dict[fr]["second_loop_length"] = second_loop.length

    potential_sequences_found_df = pd.DataFrame.from_dict(potential_sequences_found_dict, orient = "index")
    return potential_sequences_found_df


def initialize_two_loops(first_loop_dict,
                        second_loop_dict,
                        sequence_np):
    first_loop = loop()
    first_loop.initialize_from_dict(first_loop_dict)
    first_loop.fill_sequence(sequence_np)
    second_loop = loop()
    second_loop.initialize_from_dict(second_loop_dict)
    second_loop.fill_sequence(sequence_np)

    assert first_loop.check_if_complementary()
    assert second_loop.check_if_complementary()

    return first_loop, second_loop


def strengthen_one_loop(loop_to_change_orig,
                        other_loop_orig,
                        sequence_np_orig,
                        do_print = False,
                        GC_only = False):
    loop_to_change = loop_to_change_orig.copy()
    other_loop = other_loop_orig.copy()
    sequence_np = copy.deepcopy(sequence_np_orig)

    string_to_print = ""
    if do_print:
        string_to_print = ""
        string_to_print += "First loop original: \n"
        string_to_print += loop_to_change.print(do_return = True)
        string_to_print += "Second loop original: \n"
        string_to_print += other_loop.print(do_return = True)

    loop_to_change.introduce_random_sequence(GC_only = GC_only)
    loop_to_change.check_if_complementary()

    changed_sequence = loop_to_change.change_original_sequence_accordingly(sequence_np)
    other_loop.fill_sequence(changed_sequence)
    fr_unpaired = other_loop.get_fraction_unpaired(allow_canonical = True)

    if do_print:
        string_to_print += "First loop modified: \n"
        string_to_print += loop_to_change.print(do_return = True)
        string_to_print += "Second loop modified: \n"
        string_to_print += other_loop.print(do_return = True)
        string_to_print += "Fraction of unpaired nucleotides in second modified loop: %.2f" % \
                           (fr_unpaired)

    if do_print:
        print(string_to_print)

    return loop_to_change, fr_unpaired


def weaken_one_loop(loop_to_weaken_orig,
                        other_loop_orig,
                        sequence_np_orig,
                        change_left_side = True,
                        do_print = False):
    loop_to_weaken = loop_to_weaken_orig.copy()
    other_loop = other_loop_orig.copy()
    sequence_np = copy.deepcopy(sequence_np_orig)

    string_to_print = ""
    if do_print:
        string_to_print = ""
        string_to_print += "First loop original: \n"
        string_to_print += loop_to_weaken.print(do_return = True)
        string_to_print += "Second loop original: \n"
        string_to_print += other_loop.print(do_return = True)

    loop_to_weaken.change_one_side(change_left_side = change_left_side,
                                       seed=None)

    changed_sequence = loop_to_weaken.change_original_sequence_accordingly(sequence_np)
    other_loop.fill_sequence(changed_sequence)
    fr_unpaired = loop_to_weaken.get_fraction_unpaired(allow_canonical = True)

    if do_print:
        string_to_print += "First loop modified: \n"
        string_to_print += loop_to_weaken.print(do_return = True)
        string_to_print += "Second loop modified: \n"
        string_to_print += other_loop.print(do_return = True)
        string_to_print += "Fraction of unpaired nucleotides in second modified loop: %.2f" % \
                           (fr_unpaired)

    if do_print:
        print(string_to_print)

    return loop_to_weaken, other_loop, fr_unpaired



def strengthen_one_loop_wrapper(first_loop_dict,
                        second_loop_dict,
                        sequence_np,

                        const_seq_left_np,
                        const_seq_right_np,

                        min_fraction_unpaired_second_loop = 0.6,
                        max_stretch_base_pairs_full_seq = 3,

                        n_iterations = 1000,
                        n_loops_desired = 20,

                        do_print = False):

    first_loop, second_loop = initialize_two_loops(
                                        first_loop_dict,
                                        second_loop_dict,
                                        sequence_np)

    enough_loops_match = False
    matching_loops = []
    i = -1
    while (i <= n_iterations) and not enough_loops_match:
        i += 1
        # change the loop randomly
        loop_to_change, fr_unpaired = strengthen_one_loop(
                            first_loop,
                            second_loop,
                            sequence_np,
                            do_print = do_print)

        # first criteria:
        # second loop is not forming more than a fraction of its original base pairs
        if fr_unpaired < min_fraction_unpaired_second_loop:
            continue

        # get the stretches for second criteria
        max_stretch = loop_to_change.scan_original_sequence_for_bp_stretches(
            sequence_np, allow_canonical=False)
        max_stretch_left = loop_to_change.scan_other_sequence_for_bp_stretches(
            const_seq_left_np, allow_canonical=False)
        max_stretch_right = loop_to_change.scan_other_sequence_for_bp_stretches(
            const_seq_right_np, allow_canonical=False)
        max_stretch = max(max_stretch, max_stretch_left, max_stretch_right)

        # second criteria
        # the new loop can't form long paired stretches with any region of the existing sequence
        if max_stretch > max_stretch_base_pairs_full_seq:
            continue

        # a loop passed all the criteria
        matching_loops.append(loop_to_change.copy())

        if len(matching_loops) >= n_loops_desired:
            enough_loops_match = True

    return matching_loops


def weaken_one_loop_wrapper(
                        first_loop_dict,
                        second_loop_dict,
                        sequence_np,

                        const_seq_left_np,
                        const_seq_right_np,

                        min_fraction_unpaired_second_loop = 0.6,
                        max_stretch_base_pairs_full_seq = 3,

                        n_iterations = 1000,
                        n_loops_desired = 20,
                        do_print = False
                        ):
    first_loop, second_loop = initialize_two_loops(
                                        first_loop_dict,
                                        second_loop_dict,
                                        sequence_np)

    enough_loops_match = False
    matching_loops = []
    i = -1
    while (i <= n_iterations) and not enough_loops_match:
        i += 1

        # which side to change: left or right. see here: https://stackoverflow.com/questions/6824681/get-a-random-boolean-in-python/6824868
        change_left_side = random.random() >= 0.5
        #print(change_left_side)
        loop_to_weaken, other_loop, fr_unpaired = weaken_one_loop(
                        second_loop,
                        first_loop,
                        sequence_np,
                        change_left_side = change_left_side,
                        do_print=False)

        # first criteria
        # the fist loop stays unchanged
        is_other_loop_still_paired = other_loop.check_if_complementary(allow_canonical=True)
        if not is_other_loop_still_paired:
            continue

        # second criteria
        # second loop is not forming more than a fraction of its original base pairs
        if fr_unpaired < min_fraction_unpaired_second_loop:
            continue

        # get the stretches for third criteria
        max_stretch = loop_to_weaken.scan_original_sequence_for_bp_stretches_one_loop(
            sequence_np, allow_canonical = False, check_left_side = change_left_side)
        max_stretch_left = loop_to_weaken.scan_other_sequence_for_bp_stretches_one_loop(
            const_seq_left_np, allow_canonical = False, check_left_side = change_left_side)
        max_stretch_right = loop_to_weaken.scan_other_sequence_for_bp_stretches_one_loop(
            const_seq_right_np, allow_canonical = False, check_left_side = change_left_side)
        max_stretch = max(max_stretch, max_stretch_left, max_stretch_right)

        # third criteria
        # the perturbed side of the weakened loop can't form long paired stretches with any region of the existing sequence
        if max_stretch > max_stretch_base_pairs_full_seq:
            continue

        # a loop passed all the criteria
        matching_loops.append(loop_to_weaken.copy())

        if len(matching_loops) >= n_loops_desired:
            enough_loops_match = True

    return matching_loops


def list_of_loops_to_sequences(list_of_loops, sequence_np):
    strings_list = []
    for loop in list_of_loops:
        new_sequence_np = loop.change_original_sequence_accordingly(sequence_np)
        new_sequence_str = sw_utils.array_to_string(new_sequence_np)
        strings_list.append(new_sequence_str)
    return strings_list


def all_perturbations_wrapper_one_fragment(
                        conflict_dict,
                        fr_sequence,
                        const_seq_left,
                        const_seq_right,
                        min_fraction_unpaired_second_loop=0.6,
                        max_stretch_base_pairs_full_seq=3,
                        n_iterations=1000,
                        n_loops_desired=20,
                        do_print = False
                        ):
    sequence_np = sw_utils.string_to_array(fr_sequence)
    const_seq_left_np = sw_utils.string_to_array(const_seq_left)
    const_seq_right_np = sw_utils.string_to_array(const_seq_right)

    major_strengthen_list = strengthen_one_loop_wrapper(
                                    conflict_dict['major'],
                                    conflict_dict['second'],
                                    sequence_np,
                                    const_seq_left_np,
                                    const_seq_right_np,
                                    min_fraction_unpaired_second_loop = min_fraction_unpaired_second_loop,
                                    max_stretch_base_pairs_full_seq = max_stretch_base_pairs_full_seq,
                                    n_iterations = n_iterations,
                                    n_loops_desired = n_loops_desired,
                                    do_print = do_print)

    second_strengthen_list = strengthen_one_loop_wrapper(
                                    conflict_dict['second'],
                                    conflict_dict['major'],
                                    sequence_np,
                                    const_seq_left_np,
                                    const_seq_right_np,
                                    min_fraction_unpaired_second_loop = min_fraction_unpaired_second_loop,
                                    max_stretch_base_pairs_full_seq = max_stretch_base_pairs_full_seq,
                                    n_iterations = n_iterations,
                                    n_loops_desired = n_loops_desired,
                                    do_print = do_print)

    second_weaken_list = weaken_one_loop_wrapper(
                                    conflict_dict['major'],
                                    conflict_dict['second'],
                                    sequence_np,
                                    const_seq_left_np,
                                    const_seq_right_np,
                                    min_fraction_unpaired_second_loop=min_fraction_unpaired_second_loop,
                                    max_stretch_base_pairs_full_seq=max_stretch_base_pairs_full_seq,
                                    n_iterations=n_iterations,
                                    n_loops_desired=n_loops_desired,
                                    do_print=do_print)


    major_weaken_list = weaken_one_loop_wrapper(
                                    conflict_dict['second'],
                                    conflict_dict['major'],
                                    sequence_np,
                                    const_seq_left_np,
                                    const_seq_right_np,
                                    min_fraction_unpaired_second_loop=min_fraction_unpaired_second_loop,
                                    max_stretch_base_pairs_full_seq=max_stretch_base_pairs_full_seq,
                                    n_iterations=n_iterations,
                                    n_loops_desired=n_loops_desired,
                                    do_print=do_print)

    perturbations_dict = {
        "major_strengthen" : major_strengthen_list,
        "second_strengthen" : second_strengthen_list,
        "second_weaken" : second_weaken_list,
        "major_weaken" : major_weaken_list,
        "conflict" : conflict_dict,
        "sequence_np" : sequence_np
    }
    return perturbations_dict


def get_maximal_bp_stretch_that_remains_anyway(loop_to_change_orig,
                                         other_loop_orig,
                                         sequence_np_orig,
                                         do_print=False):
    # check what is the longest possible stretch that could be
    # formed by the other loop if the first loop is completely changed
    sequence_np = copy.deepcopy(sequence_np_orig)
    loop_to_change_mask = loop_to_change_orig.copy()
    other_loop_mask = other_loop_orig.copy()

    # change all nts in the first loop to None
    loop_to_change_mask.mask_out()
    masked_sequence = loop_to_change_mask.change_original_sequence_accordingly(sequence_np)
    other_loop_mask.fill_sequence(masked_sequence)
    # find the longest stretch that still remains
    longest_stretch = other_loop_mask.longest_paired_stretch(allow_canonical = True)
    # find the maximal number of base pairs that still remain
    fraction_unpaired = other_loop_mask.get_fraction_unpaired(allow_canonical=True)
    number_unparied = int(other_loop_mask.length * fraction_unpaired)
    number_paired = other_loop_mask.length - number_unparied

    return longest_stretch, number_paired


def strengthen_one_loop_criteria_2_wrapper(first_loop_dict,
                        second_loop_dict,
                        sequence_np,
                        const_seq_left_np,
                        const_seq_right_np,

                        possible_longest_remaining_stretch = 4,
                        longest_remaining_stretch_fraction = 0.5,

                        min_fraction_unpaired_second_loop = 0.6,
                        min_fraction_unpaired_second_loop_overlapping = 0.5,
                        max_stretch_base_pairs_full_seq_main = 3,
                        max_stretch_base_pairs_full_seq_backup = 4,
                        min_fraction_unpaired_full_seq_backup_1 = 0.33,
                        min_fraction_unpaired_full_seq_backup_2 = 0.25,
                        n_iterations = 100,
                        n_loops_desired = 20,
                        do_print = False):

    first_loop, second_loop = initialize_two_loops(
                                        first_loop_dict,
                                        second_loop_dict,
                                        sequence_np)

    longest_stretch, max_number_paired \
        = get_maximal_bp_stretch_that_remains_anyway(
                                first_loop,
                                second_loop,
                                sequence_np)
    max_fraction_paired_second_loop = 1 - min_fraction_unpaired_second_loop
    max_acceptable_bp_number = (max_fraction_paired_second_loop * (second_loop.length - max_number_paired)) + max_number_paired
    max_acceptable_bp_number = min(max_acceptable_bp_number, min_fraction_unpaired_second_loop_overlapping * second_loop.length)
    min_number_unpaired_full_seq_backup_1 = min_fraction_unpaired_full_seq_backup_1 * first_loop.length
    min_number_unpaired_full_seq_backup_2 = min_fraction_unpaired_full_seq_backup_2 * first_loop.length

    # first criteria: the longest possible remaining stretch should not be
    # longer than possible_longest_remaining_stretch
    if longest_stretch > possible_longest_remaining_stretch:
        return []

    # second criteria: the longest possible remaining stretch should not be
    # taking more than longest_remaining_stretch_fraction fraction of the second_loop
    if (longest_stretch / second_loop.length) > longest_remaining_stretch_fraction:
        return []

    # if the first two criteria do pass, it is possible to find a perturbation
    # so we do the search

    enough_loops_match = False
    matching_loops_main = []
    matching_loops_backups = {0 : [], 1 : [], 2 : []}
    i = -1
    while (i <= n_iterations) and not enough_loops_match:
        i += 1
        # change the loop randomly
        loop_to_change, fr_unpaired = strengthen_one_loop(
                            first_loop,
                            second_loop,
                            sequence_np,
                            do_print = do_print)

        # first & second criteria:
        # number of base pairs is not more than
        # min_fraction_unpaired_second_loop * (loop length - max_number_paired) + max_number_paired
        # and not more than min_fraction_unpaired_second_loop_overlapping
        number_unparied = int(second_loop.length * fr_unpaired)
        number_paired = second_loop.length - number_unparied
        if number_paired >  max_acceptable_bp_number:
            continue

        # get the stretches for the third criteria
        max_stretch = loop_to_change.scan_original_sequence_for_bp_stretches(
            sequence_np, allow_canonical=False)
        max_stretch_left = loop_to_change.scan_other_sequence_for_bp_stretches(
            const_seq_left_np, allow_canonical=False)
        max_stretch_right = loop_to_change.scan_other_sequence_for_bp_stretches(
            const_seq_right_np, allow_canonical=False)
        max_stretch = max(max_stretch, max_stretch_left, max_stretch_right)

        # get the number of base pairs for the third criteria
        max_number_of_base_pairs = loop_to_change.scan_original_sequence_for_number_of_bps(
            sequence_np, allow_canonical=False)
        max_number_of_base_pairs_left = loop_to_change.scan_other_sequence_for_number_of_bps(
            const_seq_left_np, allow_canonical=False)
        max_number_of_base_pairs_right = loop_to_change.scan_other_sequence_for_number_of_bps(
            const_seq_right_np, allow_canonical=False)
        max_number_of_base_pairs_overall = max(max_number_of_base_pairs, max_number_of_base_pairs_left, max_number_of_base_pairs_right)


        # third criteria
        #- the new sequence can't form paired stretches:
        # - main criteria: no paired stretches of more than 3 nt long with the other parts of the fragment
        # - backups umbrella criteria: no paired stretches of more than 4 nt long with the other parts of the fragment + the surrounding sequence
        #     - backup 1 criteria: no positions with less than 40% unpaired bases in the other parts of the fragment + the surrounding sequence
        #     - backup 2 criteria: no positions with less than 25% unpaired bases in the other parts of the fragment + the surrounding sequence
        #     - backup 3 criteria: no positions with less than 25% unpaired bases in the other parts of the fragment
        # - if no sequences match main criteria but some match the backup criteria, take those. Go down the backup criterias list until there are some fragments found

        if max_stretch <= max_stretch_base_pairs_full_seq_main:
            matching_loops_main.append(loop_to_change.copy())
        elif max_stretch <= max_stretch_base_pairs_full_seq_backup:
            min_number_of_unpaired = loop_to_change.length - max_number_of_base_pairs_overall
            min_number_of_unpaired_no_surrounding = loop_to_change.length - max_number_of_base_pairs
            if min_number_of_unpaired >= min_number_unpaired_full_seq_backup_1:
                matching_loops_backups[0].append(loop_to_change.copy())
            elif min_number_of_unpaired >= min_number_unpaired_full_seq_backup_2:
                matching_loops_backups[1].append(loop_to_change.copy())
            elif min_number_of_unpaired_no_surrounding >= min_number_unpaired_full_seq_backup_2:
                matching_loops_backups[2].append(loop_to_change.copy())


        if len(matching_loops_main) >= n_loops_desired:
            enough_loops_match = True

    if len(matching_loops_main) > 0:
        output_loops_list = matching_loops_main
    else:
        i = 0
        not_found_yet = True
        output_loops_list = []
        while (i in matching_loops_backups) and not_found_yet:
            if len(matching_loops_backups[i]) > 0:
                output_loops_list = matching_loops_backups[i]
                not_found_yet = False
            i += 1
    if len(output_loops_list) > n_loops_desired:
        output_loops_list = output_loops_list[0 : n_loops_desired]
    # print("Loop length %d, found %d, %d, %d and %d fragments, output length is " %
    #       (loop_to_change.length, len(matching_loops_main), len(matching_loops_backups[0]),
    #                                                             len(matching_loops_backups[1]),
    #                                                                 len(matching_loops_backups[2])),
    #       len(output_loops_list))

    return output_loops_list


def strengthen_one_loop_criteria_3_wrapper(first_loop_dict,
                        second_loop_dict,
                        sequence_np,
                        possible_longest_remaining_stretch = 4,
                        longest_remaining_stretch_fraction = 0.5,
                        min_fraction_unpaired_second_loop = 0.6,
                        min_fraction_unpaired_second_loop_overlapping = 0.5,
                        n_iterations = 5,
                        do_print = False):

    first_loop, second_loop = initialize_two_loops(
                                        first_loop_dict,
                                        second_loop_dict,
                                        sequence_np)

    longest_stretch, max_number_paired \
        = get_maximal_bp_stretch_that_remains_anyway(
                                first_loop,
                                second_loop,
                                sequence_np)
    max_fraction_paired_second_loop = 1 - min_fraction_unpaired_second_loop
    max_acceptable_bp_number = (max_fraction_paired_second_loop * (second_loop.length - max_number_paired)) + max_number_paired
    max_acceptable_bp_number = min(max_acceptable_bp_number, min_fraction_unpaired_second_loop_overlapping * second_loop.length)

    # first criteria: the longest possible remaining stretch should not be
    # longer than possible_longest_remaining_stretch
    if longest_stretch > possible_longest_remaining_stretch:
        return []

    # second criteria: the longest possible remaining stretch should not be
    # taking more than longest_remaining_stretch_fraction fraction of the second_loop
    if (longest_stretch / second_loop.length) > longest_remaining_stretch_fraction:
        return []

    # if the first two criteria do pass, it is possible to find a perturbation
    # so we do the search

    enough_loops_match = False
    matching_loops_main = []
    i = -1
    while i <= n_iterations:
        i += 1
        # change the loop randomly
        loop_to_change, fr_unpaired = strengthen_one_loop(
                            first_loop,
                            second_loop,
                            sequence_np,
                            GC_only = True,
                            do_print = do_print)

        # first & second criteria:
        # number of base pairs is not more than
        # min_fraction_unpaired_second_loop * (loop length - max_number_paired) + max_number_paired
        # and not more than min_fraction_unpaired_second_loop_overlapping
        number_unparied = int(second_loop.length * fr_unpaired)
        number_paired = second_loop.length - number_unparied
        if number_paired > max_acceptable_bp_number:
            continue
        matching_loops_main.append(loop_to_change.copy())
    return matching_loops_main


def weaken_one_loop_criteria_3_wrapper(
                        first_loop_dict,
                        second_loop_dict,
                        sequence_np,
                        max_fraction_unpaired_first_loop = 0.2,
                        min_fraction_unpaired_second_loop = 0.6,
                        n_iterations = 10,
                        ):
    first_loop, second_loop = initialize_two_loops(
                                        first_loop_dict,
                                        second_loop_dict,
                                        sequence_np)
    matching_loops = []
    i = -1
    while i <= n_iterations:
        i += 1

        # which side to change: left or right. see here: https://stackoverflow.com/questions/6824681/get-a-random-boolean-in-python/6824868
        change_left_side = random.random() >= 0.5
        loop_to_weaken, other_loop, fr_unpaired = weaken_one_loop(
                        second_loop,
                        first_loop,
                        sequence_np,
                        change_left_side = change_left_side,
                        do_print=False)

        # first criteria
        # the fist loop has not more than max_fraction_unpaired_changed_loop unpaired bases
        other_loop_fr_unpaired = other_loop.get_fraction_unpaired(allow_canonical=True)
        if other_loop_fr_unpaired > max_fraction_unpaired_first_loop:
            continue

        # second criteria
        # second loop is not forming more than a fraction of its original base pairs
        if fr_unpaired < min_fraction_unpaired_second_loop:
            continue

        # a loop passed all the criteria
        matching_loops.append(loop_to_weaken.copy())

    return matching_loops


def reformat_folded_perturbations_dict(folded_perturb_dict):
    out_dict = {}
    for pert_name in folded_perturb_dict:
        pert_name_split = pert_name.split('|')
        fr_name = pert_name_split[0]
        perturbation_type = pert_name_split[1]
        number = pert_name_split[2]
        curr_base_pair_array = base_pair_array(folded_perturb_dict[pert_name]["dotbracket"])
        if fr_name not in out_dict:
            out_dict[fr_name] = {}
        if perturbation_type not in out_dict[fr_name]:
            out_dict[fr_name][perturbation_type] = {}
        out_dict[fr_name][perturbation_type][number] = {}
        out_dict[fr_name][perturbation_type][number]['sequence'] = folded_perturb_dict[pert_name]['sequence']
        out_dict[fr_name][perturbation_type][number]['bp'] = curr_base_pair_array
    return out_dict


def get_perturb_combined_dictionary(potential_perturb_fn,
                                    folded_perturbations_fn):
    # upload a file with all the potential perturbations listed
    # and also a file with all the perturbations folded
    # combine these two dictionaries into one
    # note that some perturbations might be missing because of additional filtering (for example, by restriction enzyme sites
    combined_dict_loc = {}
    pert_types_list = ["major_strengthen",
                       "second_strengthen",
                       "second_weaken",
                       "major_weaken"]
    potential_perturb_dict = sw_io.read_possible_perturbation_file(potential_perturb_fn)
    folded_perturb_dict_raw = sw_io.read_folded_fasta(folded_perturbations_fn)
    folded_perturb_dict = reformat_folded_perturbations_dict(folded_perturb_dict_raw)

    for fr_name in potential_perturb_dict:
        combined_dict_loc[fr_name] = {
            "conflict": potential_perturb_dict[fr_name]['conflict'],
            "sequence_np": potential_perturb_dict[fr_name]['sequence_np'],
            "dotbracket_major": potential_perturb_dict[fr_name]['dotbracket_major'],
            "dotbracket_second": potential_perturb_dict[fr_name]['dotbracket_second']
        }
        if fr_name not in folded_perturb_dict:
            for pert_type in pert_types_list:
                combined_dict_loc[fr_name][pert_type] = {}
        else:
            for pert_type in pert_types_list:
                if pert_type not in folded_perturb_dict[fr_name]:
                    combined_dict_loc[fr_name][pert_type] = {}
                    continue
                else:
                    combined_dict_loc[fr_name][pert_type] = folded_perturb_dict[fr_name][pert_type]
    return combined_dict_loc


def does_loop_match_folding(loop, bps,
                            min_fraction_bps_matching):
    same_bps_number = loop.compare_to_base_pair_array(bps)
    same_bps_fraction = same_bps_number / loop.length
    does_folding_match = same_bps_fraction >= min_fraction_bps_matching
    return does_folding_match, same_bps_fraction


def info_on_folded_perturbs_to_dataframe(pert_type_name,
                                         inp_df):
    # take a dataframe with information on folding of several perturbations
    # return a dictionary indicating if there are any potential and/or passing perturbations
    # and containing the best passing perturbation or None if there isn't any
    potential_label = "%s_potential" % pert_type_name
    passed_label = "%s_passed" % pert_type_name
    best_one_label = "%s_best_one" % pert_type_name
    default_dict = {
        potential_label : False,
        passed_label : False,
        best_one_label : None
    }
    is_there_a_potential_pert = inp_df.shape[0] > 0
    if not is_there_a_potential_pert:
        return default_dict
    else:
        default_dict[potential_label] = True
    positive_df = inp_df[inp_df['pass']].copy()
    is_there_passed_pert = positive_df.shape[0] > 0
    if not is_there_passed_pert:
        return default_dict
    else:
        default_dict[passed_label] = True

    positive_df_sorted = positive_df.sort_values(by=['fraction_matching', 'GC'],
                                         axis=0,
                                         ascending=False)
    default_dict[best_one_label] = positive_df_sorted.iloc[0]['sequence']
    return default_dict


def do_comparisons_on_combined_dictionary(inp_dict,
                                          min_fraction_bps_matching = 0.80):
    # iterate through folded iterations in a combined dictionary from the get_perturb_combined_dictionary function
    # for each perturbation, report
    # (1) if it matches the expected loop or if it doesn't
    # (2) if the match is 100% or if it's not
    # (3) GC content
    # return a dictionary where every type of perturbations is represented by a dataframe
    out_dict = {}

    for fr_name in inp_dict:
        out_dict[fr_name] = {}
        major_loop, second_loop = initialize_two_loops(
                            inp_dict[fr_name]['conflict']['major'],
                            inp_dict[fr_name]['conflict']['second'],
                            inp_dict[fr_name]['sequence_np'])
        curr_interm_dict = {}
        for number in inp_dict[fr_name]["major_strengthen"]:
            curr_interm_dict[number] = {}
            does_folding_match, same_bps_fraction = does_loop_match_folding(
                                    major_loop,
                                    inp_dict[fr_name]["major_strengthen"][number]['bp'],
                                    min_fraction_bps_matching = min_fraction_bps_matching)
            gc_cont = sw_utils.get_GC_content(inp_dict[fr_name]["major_strengthen"][number]['sequence'])
            curr_interm_dict[number]['pass'] = does_folding_match
            curr_interm_dict[number]['fraction_matching'] = same_bps_fraction
            curr_interm_dict[number]['GC'] = gc_cont
            curr_interm_dict[number]['sequence'] = inp_dict[fr_name]["major_strengthen"][number]['sequence']
            curr_interm_dict[number]['bp'] = inp_dict[fr_name]["major_strengthen"][number]['bp']
            out_dict[fr_name]["major_strengthen"] = pd.DataFrame.from_dict(
                                                    curr_interm_dict, orient='index')

        curr_interm_dict = {}
        for number in inp_dict[fr_name]["second_strengthen"]:
            curr_interm_dict[number] = {}
            does_folding_match, same_bps_fraction = does_loop_match_folding(
                second_loop,
                inp_dict[fr_name]["second_strengthen"][number]['bp'],
                min_fraction_bps_matching=min_fraction_bps_matching)
            gc_cont = sw_utils.get_GC_content(inp_dict[fr_name]["second_strengthen"][number]['sequence'])
            curr_interm_dict[number]['pass'] = does_folding_match
            curr_interm_dict[number]['fraction_matching'] = same_bps_fraction
            curr_interm_dict[number]['GC'] = gc_cont
            curr_interm_dict[number]['sequence'] = inp_dict[fr_name]["second_strengthen"][number]['sequence']
            curr_interm_dict[number]['bp'] = inp_dict[fr_name]["second_strengthen"][number]['bp']
            out_dict[fr_name]["second_strengthen"] = pd.DataFrame.from_dict(
                                                    curr_interm_dict, orient='index')

        curr_interm_dict = {}
        for number in inp_dict[fr_name]["second_weaken"]:
            curr_interm_dict[number] = {}
            does_folding_match, same_bps_fraction = does_loop_match_folding(
                major_loop,
                inp_dict[fr_name]["second_weaken"][number]['bp'],
                min_fraction_bps_matching=min_fraction_bps_matching)
            gc_cont = sw_utils.get_GC_content(inp_dict[fr_name]["second_weaken"][number]['sequence'])
            curr_interm_dict[number]['pass'] = does_folding_match
            curr_interm_dict[number]['fraction_matching'] = same_bps_fraction
            curr_interm_dict[number]['GC'] = gc_cont
            curr_interm_dict[number]['sequence'] = inp_dict[fr_name]["second_weaken"][number]['sequence']
            curr_interm_dict[number]['bp'] = inp_dict[fr_name]["second_weaken"][number]['bp']
            out_dict[fr_name]["second_weaken"] = pd.DataFrame.from_dict(
                                                    curr_interm_dict, orient='index')

        curr_interm_dict = {}
        for number in inp_dict[fr_name]["major_weaken"]:
            curr_interm_dict[number] = {}
            does_folding_match, same_bps_fraction = does_loop_match_folding(
                second_loop,
                inp_dict[fr_name]["major_weaken"][number]['bp'],
                min_fraction_bps_matching=min_fraction_bps_matching)
            gc_cont = sw_utils.get_GC_content(inp_dict[fr_name]["major_weaken"][number]['sequence'])
            curr_interm_dict[number]['pass'] = does_folding_match
            curr_interm_dict[number]['fraction_matching'] = same_bps_fraction
            curr_interm_dict[number]['GC'] = gc_cont
            curr_interm_dict[number]['sequence'] = inp_dict[fr_name]["major_weaken"][number]['sequence']
            curr_interm_dict[number]['bp'] = inp_dict[fr_name]["major_weaken"][number]['bp']
            out_dict[fr_name]["major_weaken"] = pd.DataFrame.from_dict(
                                                    curr_interm_dict, orient='index')

    return out_dict


def extract_stats_and_best_perturb_from_combined_dictionary(inp_dict):
    # iterate through fragments in a combined dictionary from the do_comparisons_on_combined_dictionary function
    # for each fragment, report
    # (1) is there a matching perturbation per each perturbation type
    # (2) was there a potential perturbation that got denied cause it doesn't fold the way we wanted
    # (3) what is the best perturbation based on (a) fraction of matching basepairs and (b) GC content
    # return a dataframe that contains the calculated statistics and the best perturbation sequence if there is any
    pert_types_list = ["major_strengthen",
                       "second_strengthen",
                       "second_weaken",
                       "major_weaken"]
    info_dict = {}

    for fr_name in inp_dict:
        info_dict[fr_name] = {}
        for pert_type in pert_types_list:
            if pert_type in inp_dict[fr_name]:
                curr_df = inp_dict[fr_name][pert_type]
            else:
                curr_df = pd.DataFrame()
            curr_dict = info_on_folded_perturbs_to_dataframe(pert_type, curr_df)
            for col in curr_dict:
                info_dict[fr_name][col] = curr_dict[col]

    out_df = pd.DataFrame.from_dict(
                    info_dict, orient='index')

    return out_df


def processing_folded_perturbations_wrapper(potential_perturb_fn,
                                            folded_perturbations_fn,
                                            min_fraction_bps_matching = 0.80):
    combined_pert_dict_raw = get_perturb_combined_dictionary(
                                potential_perturb_fn,
                                folded_perturbations_fn)
    combined_pert_dict_processed = do_comparisons_on_combined_dictionary(
                                combined_pert_dict_raw,
                                min_fraction_bps_matching = min_fraction_bps_matching)
    combined_pert_df = extract_stats_and_best_perturb_from_combined_dictionary(
                                combined_pert_dict_processed)
    return combined_pert_df


def add_2nd_combined_pert_df_to_1st(df_1, df_2):
    pert_types_list = ["major_strengthen",
                       "second_strengthen",
                       "second_weaken",
                       "major_weaken"]
    out_df = df_1.copy()
    indexes_to_process = set(df_1.index).intersection(set(df_2.index))
    for index in indexes_to_process:
        curr_row_1 = df_1.loc[index]
        curr_row_2 = df_2.loc[index]
        for pert_type in pert_types_list:
            pass_col = "%s_passed" % pert_type
            best_one_col = "%s_best_one" % pert_type
            # if there is a perturbation that passed already - don't proceed
            if curr_row_1[pass_col]:
                continue
            # if there isn't a new passed perturbation - don't proceed
            if not curr_row_2[pass_col]:
                continue
            # at this point there is a new perturbation that passed and no old one
            out_df.at[index, pass_col] = True
            out_df.at[index, best_one_col] = curr_row_2[best_one_col]

    return out_df


def reformat_combined_pert_df(inp_df):
    out_df = pd.DataFrame(index = inp_df.index)

    out_df["strengthen_passed"] = inp_df["major_strengthen_passed"] & inp_df["second_strengthen_passed"]
    out_df["weaken_passed"] = inp_df["second_weaken_passed"] & inp_df["major_weaken_passed"]
    out_df["too_much_overlap"] = inp_df["both_weakens_fail"] & inp_df["both_sides_overlap"]

    out_df["major_strengthen_pert"] = inp_df["major_strengthen_best_one"]
    out_df["second_strengthen_pert"] = inp_df["second_strengthen_best_one"]
    out_df["major_weaken_pert"] = inp_df["major_weaken_best_one"]
    out_df["second_weaken_pert"] = inp_df["second_weaken_best_one"]
    return out_df


def do_both_sides_overlap(first_loop, second_loop):
    left_overlap = sw_utils.get_overlap(
        (first_loop.left_start, first_loop.left_end),
        (second_loop.left_start, second_loop.left_end))
    right_overlap = sw_utils.get_overlap(
        (first_loop.right_start, first_loop.right_end),
        (second_loop.right_start, second_loop.right_end))

    return (left_overlap > 0) & (right_overlap > 0)


def check_if_both_sides_overlap_on_perturb_dict_item(dict_item):
    first_loop, second_loop = initialize_two_loops(
                                dict_item['conflict']['major'],
                                dict_item['conflict']['second'],
                                dict_item["sequence_np"])
    do_both_overlap = do_both_sides_overlap(first_loop, second_loop)
    return do_both_overlap


def get_mass_centers_for_mut_combination(source_df,
                                   df_with_values,
                                   pert_types_list=[
                                                  "orig",
                                                  "major_strengthen_pert",
                                                  "second_strengthen_pert",
                                                  "second_weaken_pert",
                                                  "major_weaken_pert"]
                                         ):
    mass_center_values_df_1 = pd.DataFrame(index = source_df.index,
                                         columns = pert_types_list)
    mass_center_values_df_2 = pd.DataFrame(index = source_df.index,
                                         columns = pert_types_list)

    for index, source_row in source_df.iterrows():
        general_sub_df = df_with_values[(df_with_values["fr_name"] == source_row["fr_name"]) &
                                        # or original
                                        ((df_with_values["source"] == source_row["source"]) |
                                         (df_with_values["source"] == "orig"))]
        for pert_type in pert_types_list:
            sub_df = general_sub_df[general_sub_df['pert_type'] == pert_type]
            if sub_df.shape[0] == 0:
                continue
            assert sub_df.shape[0] == 1
            value_1 = sub_df['mass_center_1'].values.tolist()[0]
            value_2 = sub_df['mass_center_2'].values.tolist()[0]
            mass_center_values_df_1.at[index, pert_type] = value_1
            mass_center_values_df_2.at[index, pert_type] = value_2

    return mass_center_values_df_1, mass_center_values_df_2


def mutations_5_differ_strict_rule(row):
    reg_phenotype = False
    if row['major_strengthen_pert'] > row['second_strengthen_pert'] and \
        row['major_strengthen_pert'] > row['major_weaken_pert'] and \
        row['second_weaken_pert'] > row['second_strengthen_pert'] and \
        row['second_weaken_pert'] > row['major_weaken_pert']:
        # major conformation is more active than second conformation
        reg_phenotype = True
    if row['second_strengthen_pert'] > row['major_strengthen_pert'] and \
        row['second_strengthen_pert'] > row['second_weaken_pert'] and \
        row['major_weaken_pert'] > row['major_strengthen_pert'] and \
        row['major_weaken_pert'] > row['second_weaken_pert']:
        # second conformation is more active than major conformation
        reg_phenotype = True
    return reg_phenotype


def mutations_5_differ_orig_strict_rule(row):
    reg_phenotype = False
    if row['major_strengthen_pert'] > row['second_strengthen_pert'] and \
        row['major_strengthen_pert'] > row['major_weaken_pert'] and \
        row['second_weaken_pert'] > row['second_strengthen_pert'] and \
        row['second_weaken_pert'] > row['major_weaken_pert'] and \
        row['orig'] < row['major_strengthen_pert'] and \
        row['orig'] < row['second_weaken_pert'] and \
        row['orig'] > row['major_weaken_pert'] and \
        row['orig'] > row['second_strengthen_pert']:
        # major conformation is more active than second conformation
        # original fragment is in between
        reg_phenotype = True
    if row['second_strengthen_pert'] > row['major_strengthen_pert'] and \
        row['second_strengthen_pert'] > row['second_weaken_pert'] and \
        row['major_weaken_pert'] > row['major_strengthen_pert'] and \
        row['major_weaken_pert'] > row['second_weaken_pert'] and \
        row['orig'] < row['major_weaken_pert'] and \
        row['orig'] < row['second_strengthen_pert'] and \
        row['orig'] > row['major_strengthen_pert'] and \
        row['orig'] > row['second_weaken_pert']:
        # second conformation is more active than major conformation
        # original fragment is in between
        reg_phenotype = True
    return reg_phenotype


def subset_counts_five_perturbations(source_row,
                                     df_with_values,
                                     pert_types_order,
                                     replicate_1_column_names,
                                     replicate_2_column_names,
                                     return_subset_order = False
                                     ):
    general_sub_df = df_with_values[(df_with_values["fr_name"] == source_row["fr_name"]) &
                                    # or original
                                    ((df_with_values["source"] == source_row["source"]) |
                                     (df_with_values["source"] == "orig"))].copy()
    general_sub_df.index = general_sub_df['pert_type']
    general_sub_df_ordered = general_sub_df.loc[pert_types_order]
    rep_1_counts = general_sub_df_ordered[replicate_1_column_names]
    rep_2_counts = general_sub_df_ordered[replicate_2_column_names]
    return rep_1_counts, rep_2_counts


def subset_counts_five_perturbations_from_averaged(source_row,
                                     df_with_values,
                                     pert_types_order,
                                     column_names,
                                     ):
    general_sub_df = df_with_values[(df_with_values["fr_name"] == source_row["fr_name"]) &
                                    # or original
                                    ((df_with_values["source"] == source_row["source"]) |
                                     (df_with_values["source"] == "orig"))].copy()
    general_sub_df.index = general_sub_df['pert_type']
    general_sub_df_ordered = general_sub_df.loc[pert_types_order]
    counts = general_sub_df_ordered[column_names]
    return counts


def calculate_correlations_for_perturbations(pert_count_df,
                                             method = 'pearson'):
    return pert_count_df.transpose().corr(method = method)


def average_correlation_for_similar_perturbations_5(pert_corr_df,
                                                    geom_mean = True):
    sim_exp_1 = pert_corr_df.loc['major_strengthen_pert', 'second_weaken_pert']
    sim_exp_2 = pert_corr_df.loc['second_strengthen_pert', 'major_weaken_pert']
    dissim_exp_1 = pert_corr_df.loc['major_strengthen_pert', 'second_strengthen_pert']
    dissim_exp_2 = pert_corr_df.loc['major_strengthen_pert', 'major_weaken_pert']
    dissim_exp_3 = pert_corr_df.loc['second_weaken_pert', 'second_strengthen_pert']
    dissim_exp_4 = pert_corr_df.loc['second_weaken_pert', 'major_weaken_pert']
    if geom_mean:
        # ig using geometric mean, add 2 to make all the correlation values positive
        sim_avg = sw_utils.geom_mean_1d(np.array([
                                        sim_exp_1 + 2,
                                        sim_exp_2 + 2]))
        dissim_avg = sw_utils.geom_mean_1d(np.array([
                                        dissim_exp_1 + 2,
                                        dissim_exp_2 + 2,
                                        dissim_exp_3 + 2,
                                        dissim_exp_4 + 2]))
    else:
        sim_avg = (sim_exp_1 + sim_exp_2) / 2
        dissim_avg = (dissim_exp_1 + dissim_exp_2 + dissim_exp_3 + dissim_exp_4) / 4
    return sim_avg, dissim_avg


def choose_fragments_by_indices(inp_df, indices_set):
    boolean_mask = inp_df.apply(lambda x:
                                True if "%s_%s" % (x['fr_name'], x['source']) in indices_set else False,
                                axis = 1)
    out_df = inp_df[boolean_mask].copy()
    return out_df


def separate_fragments_given_number_of_mutations(inp_df,
                                                 n_mutations,
                                                 do_print = True):
    inp_df_perturb_counts_df = inp_df.groupby(['fr_name', 'source']).size().reset_index().\
                            rename({0:"count"}, axis = 1)
    out_df = inp_df_perturb_counts_df[inp_df_perturb_counts_df['count'] == n_mutations]
    out_df.index = out_df.apply(lambda x: "%s|%s" % (x['fr_name'], x['source']), axis = 1)
    if do_print:
        print("Number of fragments with %d mutations is %d" % (n_mutations, out_df.shape[0]))
    return out_df

def make_mutations_summary(inp_df):
    inp_df_perturb_counts_df = inp_df.groupby(['fr_name', 'source']).size().reset_index().\
                            rename({0:"count"}, axis = 1)
    out_df = inp_df_perturb_counts_df[inp_df_perturb_counts_df['count'] > 1]
    out_df.index = out_df.apply(lambda x: "%s|%s" % (x['fr_name'], x['source']), axis = 1)
    return out_df

def calculate_correlations_for_5_perturbations(mutations_df,
                                               counts_df,
                                               perturbations_order,
                                               replicate_1_column_names,
                                               replicate_2_column_names,
                                               geom_mean = True
                                               ):
    corr_df = pd.DataFrame(index=mutations_df.index,
                           columns=['fr_name', 'source',
                                    'similar_avg_1', 'dissimilar_avg_1',
                                    'similar_avg_2', 'dissimilar_avg_2',
                                    "similar_avg_both",
                                    "dissimilar_avg_both",
                                    "avg_difference_1",
                                    "avg_difference_2",
                                    "avg_difference_both"
                                    ])

    corr_df["fr_name"] = mutations_df['fr_name']
    corr_df["source"] = mutations_df['source']
    corr_df["gene_name"] = mutations_df['gene_name']

    for index, row in mutations_df.iterrows():
        counts_rep_1, counts_rep_2 = subset_counts_five_perturbations(
            row,
            counts_df,
            perturbations_order,
            replicate_1_column_names,
            replicate_2_column_names)
        corr_mtx_1 = calculate_correlations_for_perturbations(counts_rep_1)
        corr_mtx_2 = calculate_correlations_for_perturbations(counts_rep_2)
        sim_avg_1, dissim_avg_1 = average_correlation_for_similar_perturbations_5(corr_mtx_1, geom_mean)
        sim_avg_2, dissim_avg_2 = average_correlation_for_similar_perturbations_5(corr_mtx_2, geom_mean)
        corr_df.at[index, 'similar_avg_1'] = sim_avg_1
        corr_df.at[index, 'similar_avg_2'] = sim_avg_2
        corr_df.at[index, 'dissimilar_avg_1'] = dissim_avg_1
        corr_df.at[index, 'dissimilar_avg_2'] = dissim_avg_2
        if geom_mean:
            corr_df.at[index, "similar_avg_both"] = sw_utils.geom_mean_1d(np.array([
                                                    sim_avg_1,
                                                    sim_avg_2]))
            corr_df.at[index, "dissimilar_avg_both"] = sw_utils.geom_mean_1d(np.array([
                                                    dissim_avg_1,
                                                    dissim_avg_2]))
        else:
            corr_df.at[index, "similar_avg_both"] = (sim_avg_1 + sim_avg_2) / 2
            corr_df.at[index, "dissimilar_avg_both"] = (dissim_avg_1 + dissim_avg_2) / 2
        corr_df.at[index, "avg_difference_1"] = sim_avg_1 - dissim_avg_1
        corr_df.at[index, "avg_difference_2"] = sim_avg_2 - dissim_avg_2
        if geom_mean:
            corr_df.at[index, "avg_difference_both"] = sw_utils.geom_mean_1d(np.array([
                                                    corr_df.loc[index, "avg_difference_1"] + 2,
                                                    corr_df.loc[index, "avg_difference_2"] + 2]))
        else:
            corr_df.at[index, "avg_difference_both"] = ((sim_avg_1 - dissim_avg_1) + (sim_avg_2 - dissim_avg_2)) / 2

    return corr_df


def print_counts_one_fragment(inp_df, fr_name, include_barcode = False):
    columns_of_interest = ['fr_name', 'source', 'pert_type']
    if include_barcode:
        columns_of_interest.append('barcode')
    columns_of_interest += ['1_1', '1_2', '1_3', '1_4', '1_5', '1_6', '1_7', '1_8',
        '2_1', '2_2', '2_3', '2_4', '2_5', '2_6', '2_7', '2_8']
    sub_df = inp_df[inp_df['fr_name'] == fr_name][columns_of_interest]
    return sub_df


def fragment_name_to_gene_names(transcript_id,
                        transcripts_to_genes_dict,
                        ensembl_to_gene_names_dict,
                        first_dict_direct = False):
    transcript_id_short = transcript_id.split('_')[0]
    if transcript_id_short not in transcripts_to_genes_dict:
        return ''
    gene_id = transcripts_to_genes_dict[transcript_id_short]
    if first_dict_direct:
        return gene_id
    if gene_id not in ensembl_to_gene_names_dict:
        return ''
    gene_name = ensembl_to_gene_names_dict[gene_id]
    return gene_name


def add_gene_names_to_df(inp_df,
                         transcripts_to_genes_dict,
                         ensembl_to_gene_names_dict,
                         first_dict_direct = False):
    out_df = inp_df.copy()
    out_df['gene_name'] = out_df.apply(lambda x:
                                       fragment_name_to_gene_names(x['fr_name'],
                                                                   transcripts_to_genes_dict,
                                                                   ensembl_to_gene_names_dict,
                                                                   first_dict_direct = first_dict_direct),
                                       axis = 1
                                       )
    return out_df




def fit_regression_one_fragment_one_perturbations(
                                            row,
                                            dna_counts_df,
                                            rna_counts_df,
                                            perturbations_order,
                                            replicate_1_column_names,
                                            replicate_2_column_names,
                                            leave_out = 'orig',
                                            n_bins = 8,
                                            n_reps = 2,
                                            do_print_formula = False
                                            ):
    # get counts arrays
    dna_counts_rep_1, dna_counts_rep_2 = subset_counts_five_perturbations(
                                        row,
                                        dna_counts_df,
                                        perturbations_order,
                                        replicate_1_column_names,
                                        replicate_2_column_names)
    rna_counts_rep_1, rna_counts_rep_2 = subset_counts_five_perturbations(
                                        row,
                                        rna_counts_df,
                                        perturbations_order,
                                        replicate_1_column_names,
                                        replicate_2_column_names)
    # merge the replicate measurements together
    dna_counts_concat = np.concatenate((dna_counts_rep_1, dna_counts_rep_2), axis=1)
    rna_counts_concat = np.concatenate((rna_counts_rep_1, rna_counts_rep_2), axis=1)
    # flatten the counts arrays
    dna_counts_flatten = dna_counts_concat.flatten()
    rna_counts_flatten = rna_counts_concat.flatten()
    # merge DNA and RNA
    values_array = np.concatenate((dna_counts_flatten, rna_counts_flatten))

    # make a dataframe with covariates
    # 4 perturbations + original fragment
    # describing N = 5 classes with N - 1 = 4 covariates
    # for each perturbation, I have 16 measurements for DNA and 16 measurements for RNA
    n_perturbations = len(perturbations_order)
    n_samples = n_bins * n_reps
    total_measurements = n_samples * n_perturbations
    covariate_template = np.zeros(total_measurements, dtype = bool)
    covariates_array = np.zeros((total_measurements * 2, n_perturbations),
                                dtype = bool)
    for i, pert_name in enumerate(perturbations_order):
        curr_covariate = covariate_template.copy()
        beginning = i * n_samples
        end = (i + 1) * n_samples
        curr_covariate[beginning : end] = True
        # double the covariate since we are putting both DNA and RNA counts in
        curr_covariate_doubled = np.concatenate((curr_covariate, curr_covariate))
        covariates_array[ : , i] = curr_covariate_doubled
    covariates_df = pd.DataFrame(covariates_array,
                                 columns = perturbations_order)

    # make response variable
    # I have 5 perturbations and 16 measurements for each of them,
    # so 80 samples for RNA and 80 samples for DNA
    # first goes all the RNA, then all the DNA
    response_variable = np.array([0, 1], dtype=np.float)
    response_variable = np.repeat(response_variable, (total_measurements, total_measurements))
    intersept = np.ones_like(response_variable)
    covariates_df['intersept'] = intersept
    covariates_df['response'] = response_variable

    # make the formula
    allowed_perturbations_list = [x for x in perturbations_order if x != leave_out]
    formula_covariates = " + ".join(allowed_perturbations_list)
    formula = "%s ~ %s + 1" % ('response', formula_covariates)
    if do_print_formula:
        print(formula)

    glm_fit = GLM.from_formula(
                   formula = formula,
                   data = covariates_df,
                   family = sm.families.Binomial(link = sm.genmod.families.links.logit),
                   freq_weights = values_array
                 ).fit()

    return glm_fit


def contrast_analysis_stability_glm(glm_fit):
    # based on these two answers
    # https://stackoverflow.com/questions/34231016/compare-contrasts-in-linear-model-in-python-like-rs-contrast-library
    # https://www.statsmodels.org/dev/generated/statsmodels.discrete.discrete_model.LogitResults.t_test.html
    sim_couple_1 = ('major_strengthen_pert', 'second_weaken_pert')
    sim_couple_2 = ('second_strengthen_pert', 'major_weaken_pert')
    dissim_couple_1 = ('major_strengthen_pert', 'second_strengthen_pert')
    dissim_couple_2 = ('major_strengthen_pert', 'major_weaken_pert')
    dissim_couple_3 = ('second_weaken_pert', 'second_strengthen_pert')
    dissim_couple_4 = ('second_weaken_pert', 'major_weaken_pert')

    couples_list = [sim_couple_1, sim_couple_2, dissim_couple_1, dissim_couple_2, dissim_couple_3, dissim_couple_4]
    couples_names_list = ["sim_couple_1", "sim_couple_2", "dissim_couple_1", "dissim_couple_2", "dissim_couple_3",
                          "dissim_couple_4"]
    comparisons_list = ["%s[T.True] - %s[T.True]" % (x[0], x[1]) for x in couples_list]
    t_test = glm_fit.t_test(comparisons_list)
    summary = t_test.summary().as_csv().replace(' ','')
    summary_df = pd.read_csv(io.StringIO(summary), sep = ',', skiprows = 1, index_col = 0)
    summary_df.index = couples_names_list

    return summary_df


def geometric_means_contrast_zscores(contrast_table):
    contrast_table_copy = contrast_table.copy()
    contrast_table_copy['difference_size'] = contrast_table_copy['z'].abs()
    similar_z_scores = contrast_table_copy.loc[
                            ["sim_couple_1", "sim_couple_2"],
                            'difference_size'].to_numpy()
    dissimilar_z_scores = contrast_table_copy.loc[
                            ['dissim_couple_1', 'dissim_couple_2',
                             'dissim_couple_3', 'dissim_couple_4'],
                            'difference_size'].to_numpy()
    similar_z_scores_geom_mean = sw_utils.geom_mean_1d(similar_z_scores)
    dissimilar_z_scores_geom_mean = sw_utils.geom_mean_1d(dissimilar_z_scores)
    return similar_z_scores_geom_mean, dissimilar_z_scores_geom_mean


def log_regression_all_perturbations(
                                   mutations_df, # so that i iterate through original fragments and not perturbations
                                   gDNA_counts_df,
                                   RNA_counts_df,
                                   perturbations_order,
                                   replicate_1_column_names,
                                   replicate_2_column_names,
                                   do_print = True
                                   ):
    tic = time.time()
    stability_df = pd.DataFrame(index=mutations_df.index,
                           columns=['fr_name', 'source', 'gene_name',
                                    'similar_score', 'dissimilar_score', 'difference'])
    stability_df["fr_name"] = mutations_df['fr_name']
    stability_df["source"] = mutations_df['source']
    stability_df["gene_name"] = mutations_df['gene_name']

    for index, row in mutations_df.iterrows():
        glm_fit = fit_regression_one_fragment_one_perturbations(
            row,
            gDNA_counts_df,
            RNA_counts_df,
            perturbations_order,
            replicate_1_column_names,
            replicate_2_column_names
        )
        contrast_table = contrast_analysis_stability_glm(glm_fit)
        sim_score, dissim_score = geometric_means_contrast_zscores(contrast_table)
        stability_df.at[index, 'similar_score'] = sim_score
        stability_df.at[index, 'dissimilar_score'] = dissim_score
        stability_df.at[index, 'difference'] = dissim_score - sim_score

    stability_df = stability_df.sort_values(by='difference', ascending=False)


    toc = time.time()
    if do_print:
        print("Regression done, it took %d seconds" % (toc-tic))

    return stability_df


def transform_selection_dataframe(
                             inp_filename,
                             values_dictionary):
    curr_df = sw_io.read_selection_dataframe(inp_filename)
    curr_df.index = curr_df.apply(lambda x: "%s|%s" % (x['fragment'], x['source']), axis = 1)
    for column in curr_df.columns.tolist():
        curr_df[column] = curr_df.apply(lambda x:
                                        x[column] if x[column] not in values_dictionary \
                                            else values_dictionary[x[column]],
                                        axis = 1)
    return curr_df


def list_all_perturbations_given_loop(
                        seq,
                        loop,
                        length = 2,
                        min_from_end = 1,
                        change_left_side = True,
                        name_prefix = ""
                        ):
    mutated_seq_tuples = []

    # prepare the necessary variables
    sequence_np = sw_utils.string_to_array(seq)
    if change_left_side:
        start = loop.left_start
        end = loop.left_end
    else:
        start = loop.right_start
        end = loop.right_end

    all_combination_tuples_list = list(itertools.product(glob_vars.nt_list, repeat = length))
    all_combinations_list = [np.array(x, dtype = np.uint8) for x in all_combination_tuples_list]

    for i in range((start + min_from_end), (end - min_from_end - length + 1)):
        for subseq in all_combinations_list:
            # make sure all the nucleotides got replaced; otherwise skip
            if not (sequence_np[i : i + length] != subseq).all():
                continue
            rescue_loop = loop.copy()
            new_seq_array = sequence_np.copy()
            new_seq_array[i : i + length] = subseq
            rescue_loop.change_one_side_given_seq(
                                  new_seq_array[start : end],
                                  left_subset_border = i - start,
                                  right_subset_border = i - start + length,
                                  change_left_side = change_left_side,
                                  make_complementary = True)
            rescue_seq_array = rescue_loop.change_original_sequence_accordingly(new_seq_array)
            new_seq = sw_utils.array_to_string(new_seq_array)
            rescue_seq = sw_utils.array_to_string(rescue_seq_array)
            mutation_name = "%s_%d_%s" % (name_prefix, i, sw_utils.array_to_string(subseq))
            mutated_seq_tuples.append((mutation_name, new_seq, rescue_seq))
    return mutated_seq_tuples


def mutated_tuples_into_its_own_df(
                            inp_df,
                            curr_column_name,
                            ):
    combined_df = pd.DataFrame(columns =
                               ["gene_name", "fr_name", "source",
                                "first_loop", "second_loop", "changed_loop",
                                "mut_name", "original", "mutated", "rescued"])

    for index, row in inp_df.iterrows():
        mutated_seq_tuples = inp_df.loc[index, curr_column_name]
        for tup in mutated_seq_tuples:
            name, curr_mut, curr_resc = tup
            combined_df = combined_df.append({
                                        "gene_name" : inp_df.loc[index, 'gene_name'],
                                        "fr_name" : inp_df.loc[index, 'fr_name'],
                                        "source" : inp_df.loc[index, 'source'],
                                        "first_loop" : inp_df.loc[index, 'first_loop'],
                                        "second_loop" : inp_df.loc[index, 'second_loop'],
                                        # what loop is being changed should be written in the name
                                        "changed_loop" : name.split('_')[0],
                                        "mut_name" : name,
                                        "original" : inp_df.loc[index, 'seq'],
                                        "mutated" : curr_mut,
                                        "rescued" : curr_resc},
                                    ignore_index = True)

    return combined_df


def gather_mutations_into_a_single_dataframe(selection_df,
                                             mut_names_array):
    # for each mutation, I will keep:
    # - gene name
    # - fragment name and source (for fragment identification)
    # - the original sequence, both mutated and rescued sequence
    # - both loops (so that I know their coordinates when checking for correct folding)
    # - the name of the mutation of interest
    combined_df = pd.DataFrame(columns =
                               ["gene_name", "fr_name", "source",
                                "first_loop", "second_loop", "changed_loop",
                                "mut_name", "original", "mutated", "rescued"])

    for mut_name in mut_names_array:
        curr_df = mutated_tuples_into_its_own_df(
                            selection_df,
                            mut_name)
        combined_df = combined_df.append(curr_df, ignore_index=True)

    return combined_df


def get_mutation_bp_probabilities(
                    fr_name,
                    fr_source,
                    mibp_dict,
                    shape_dict,
                    sequence,
                    mut_sequence,
                    rescued_sequence,
                    temp_files_folder,
                    RNAstructure_path
                    ):
    matrices_constrained_tuple = folding_api.launch_folding_get_probabilities_RNAstructure(
                        fr_name,
                        sequence,
                        mibp_dict[fr_source][fr_name],
                        shape_array = shape_dict[fr_source][fr_name],
                        temp_files_folder = temp_files_folder,
                        no_constraints = False,
                        RNAstructure_path = RNAstructure_path
                        )
    matrices_original = folding_api.launch_folding_get_probabilities_RNAstructure(
                        fr_name,
                        sequence,
                        mibp_dict[fr_source][fr_name],
                        shape_array = shape_dict[fr_source][fr_name],
                        temp_files_folder = temp_files_folder,
                        no_constraints = True,
                        RNAstructure_path = RNAstructure_path
                        )
    matrix_mutated = folding_api.launch_folding_get_probabilities_RNAstructure(
                        fr_name,
                        mut_sequence,
                        mibp_dict[fr_source][fr_name],
                        shape_array = shape_dict[fr_source][fr_name],
                        temp_files_folder = temp_files_folder,
                        no_constraints = True,
                        RNAstructure_path = RNAstructure_path
                        )
    matrix_rescued = folding_api.launch_folding_get_probabilities_RNAstructure(
                        fr_name,
                        rescued_sequence,
                        mibp_dict[fr_source][fr_name],
                        shape_array = None,
                        temp_files_folder = temp_files_folder,
                        no_constraints = True,
                        RNAstructure_path = RNAstructure_path
                        )
    out_matrix_tuple = (matrices_constrained_tuple[0],
                        matrices_constrained_tuple[1],
                        matrices_original,
                        matrix_mutated,
                        matrix_rescued)
    return out_matrix_tuple


def get_probabilities_apply_to_mutation_df(
                                    mutations_df,
                                    mibp_dict,
                                    shape_dict,
                                    temp_files_folder,
                                    RNAstructure_path,
                                    n_processes = 2
                                    ):
    inp_df = mutations_df.copy()
    inp_df['matrix_tuple'] = inp_df.apply(lambda x:
                                          get_mutation_bp_probabilities(
                                              x['fr_name'],
                                              x['source'],
                                              mibp_dict,
                                              shape_dict,
                                              x['original'],
                                              x['mutated'],
                                              x['rescued'],
                                              temp_files_folder,
                                              RNAstructure_path
                                          ),
                                          axis = 1)
    out_df = mutations_df.copy()
    out_df["matrix_constrained_major"] = inp_df.apply(lambda x: x['matrix_tuple'][0], axis = 1)
    out_df["matrix_constrained_second"] = inp_df.apply(lambda x: x['matrix_tuple'][1], axis = 1)
    out_df["matrix_original"] = inp_df.apply(lambda x: x['matrix_tuple'][2], axis = 1)
    out_df["matrix_mutated"] = inp_df.apply(lambda x: x['matrix_tuple'][3], axis = 1)
    out_df["matrix_rescued"] = inp_df.apply(lambda x: x['matrix_tuple'][4], axis = 1)

    return out_df


def get_loop(source, fr_name, loop_name, sequence, MIBP_general_dict):
    first_loop = loop()
    first_loop.initialize_from_dict(MIBP_general_dict[source][fr_name][loop_name])
    sequence_np = sw_utils.string_to_array(sequence)
    first_loop.fill_sequence(sequence_np)
    return first_loop


def get_loop_exceptions(source, fr_name, loop_name, sequence, MIBP_general_dict):
    first_loop = loop()
    # for exceptions like ENST00000263388_part_9_avg_original
    if "_avg" in fr_name:
        fr_name = fr_name.replace("_avg", "")
    elif "_rep_1" in fr_name:
        fr_name = fr_name.replace("_rep_1", "")
    elif "_rep_2" in fr_name:
        fr_name = fr_name.replace("_rep_2", "")
    first_loop.initialize_from_dict(MIBP_general_dict[source][fr_name][loop_name])
    sequence_np = sw_utils.string_to_array(sequence)
    first_loop.fill_sequence(sequence_np)
    return first_loop


def extract_loop_probabilities_from_bp_matrix(loop, bp_matrix):
    prababilities_array = np.zeros(loop.length)
    for i, left, right in zip(range(loop.length),
                              range(loop.left_start, loop.left_end),
                              range(loop.right_end - 1, loop.right_start - 1, -1)):
        prababilities_array[i] = bp_matrix[left, right]

    return prababilities_array


def extract_average_loop_probability_from_bp_matrix(loop, bp_matrix, how = "mean"):
    prababilities_array = extract_loop_probabilities_from_bp_matrix(loop, bp_matrix)
    if how == "mean":
        average_prob = np.mean(prababilities_array)
    elif how == "median":
        average_prob = np.median(prababilities_array)
    else:
        sys.exit("Unknown keyword!")
    return average_prob


def extract_all_average_loop_probabilities_single_mutation(matrix_df,
                                                           how = "mean"):
    out_df = matrix_df.copy()
    for matrix_name in ["matrix_constrained_major",
                        "matrix_constrained_second",
                        "matrix_original",
                        "matrix_mutated",
                        "matrix_rescued"]:
        for loop_name in ['first_loop', 'second_loop']:
            score_name = "%s_%s" % (matrix_name.replace("matrix_", ""), loop_name)
            out_df[score_name] = out_df.apply(lambda x:
                                      extract_average_loop_probability_from_bp_matrix(x[loop_name],
                                                                                x[matrix_name],
                                                                                how = how),
                                              axis = 1)
    return out_df


def mutation_changed_bp_ratios(matrix_prob_df):
    out_df = matrix_prob_df.copy()
    out_df['mutated_first_loop_ratio'] = out_df["mutated_first_loop"] / out_df['original_first_loop']
    out_df['mutated_second_loop_ratio'] = out_df["mutated_second_loop"] / out_df['original_second_loop']
    out_df['rescued_first_loop_ratio'] = out_df["rescued_first_loop"] / out_df['original_first_loop']
    out_df['rescued_second_loop_ratio'] = out_df["rescued_second_loop"] / out_df['original_second_loop']
    return out_df


def filter_non_disruptive_mutations(mut_prob_df,
                                    max_probability_disrupted = 0.1,
                                    do_allow_coefficient = True,
                                    coefficient = 0.5,
                                    fraction_filtered_for_printing=0.8,
                                    do_print = True):
    filtered_df = pd.DataFrame(columns = mut_prob_df.columns)
    for i, group in mut_prob_df.groupby(["gene_name", "changed_loop"]):
        changed_loop = i[1]
        if changed_loop == 'major':
            original_loop_probability = group['original_first_loop'].tolist()[0]
            # probability of the mutated loop shouldn't be more than the probability of nonmutated loop
            # multiplied by some coefficient
            if do_allow_coefficient:
                max_allowed_probability = min(max_probability_disrupted, original_loop_probability * coefficient)
            else:
                max_allowed_probability = max_probability_disrupted
            probability_low_enough = group['mutated_first_loop'] <= max_allowed_probability
        elif changed_loop == 'second':
            original_loop_probability = group['original_second_loop'].tolist()[0]
            # probability of the mutated loop shouldn't be more than the probability of nonmutated loop
            # multiplied by some coefficient
            if do_allow_coefficient:
                max_allowed_probability = min(max_probability_disrupted, original_loop_probability * coefficient)
            else:
                max_allowed_probability = max_probability_disrupted
            probability_low_enough = group['mutated_second_loop'] <= max_allowed_probability
        else:
            sys.exit("unknown changed loop!")

        fraction_filtered = probability_low_enough.sum() / probability_low_enough.shape[0]
        if fraction_filtered < fraction_filtered_for_printing:
            if do_print:
                print("%s %s: kept %.2f sequences (%d out of %d total)" %
                  (i[0], i[1], fraction_filtered, probability_low_enough.sum(), probability_low_enough.shape[0]))

        filtered_df = pd.concat([filtered_df, group[probability_low_enough]])
    if do_print:
        print("For all the other ones, kept at least %.2f of all the sequences" % fraction_filtered_for_printing)

    return filtered_df


def filter_bad_rescue_mutations(mut_prob_df,
                                min_probablity_rescued = 0.4,
                                do_allow_coefficient=True,
                                coefficient = 0.8,
                                fraction_filtered_for_printing = 0.8,
                                    do_print = True):
    filtered_df = pd.DataFrame(columns = mut_prob_df.columns)
    for i, group in mut_prob_df.groupby(["gene_name", "changed_loop"]):
        changed_loop = i[1]
        if changed_loop == 'major':
            original_loop_probability = group['original_first_loop'].tolist()[0]
            # probability of the mutated loop shouldn't be more than the probability of nonmutated loop
            # multiplied by some coefficient
            if do_allow_coefficient:
                min_allowed_probability = min(min_probablity_rescued, original_loop_probability * coefficient)
            else:
                min_allowed_probability = min_probablity_rescued
            probability_high_enough = group['rescued_first_loop'] >= min_allowed_probability
        elif changed_loop == 'second':
            original_loop_probability = group['original_second_loop'].tolist()[0]
            # probability of the mutated loop shouldn't be more than the probability of nonmutated loop
            # multiplied by some coefficient
            if do_allow_coefficient:
                min_allowed_probability = min(min_probablity_rescued, original_loop_probability * coefficient)
            else:
                min_allowed_probability = min_probablity_rescued
            probability_high_enough = group['rescued_second_loop'] >= min_allowed_probability
        else:
            sys.exit("unknown changed loop!")

        fraction_filtered = probability_high_enough.sum() / probability_high_enough.shape[0]
        if fraction_filtered < fraction_filtered_for_printing:
            if do_print:
                print("%s %s: kept %.2f sequences (%d out of %d total)" %
                  (i[0], i[1], fraction_filtered, probability_high_enough.sum(), probability_high_enough.shape[0]))

        filtered_df = pd.concat([filtered_df, group[probability_high_enough]])

    return filtered_df


def filter_double_disruptive_mutations(mut_prob_df,
                                       negligible_probability = 0.1,
                                       fraction_filtered_for_printing = 0.8,
                                    do_print = True):
    filtered_df = pd.DataFrame(columns = mut_prob_df.columns)
    for i, group in mut_prob_df.groupby(["gene_name", "changed_loop"]):
        either_loop_present = (group['mutated_first_loop'] > negligible_probability) | \
                                 (group['mutated_second_loop'] > negligible_probability)

        fraction_filtered = either_loop_present.sum() / either_loop_present.shape[0]
        if fraction_filtered < fraction_filtered_for_printing:
            if do_print:
                print("%s %s: kept %.2f sequences (%d out of %d total)" %
                  (i[0], i[1], fraction_filtered, either_loop_present.sum(), either_loop_present.shape[0]))

        filtered_df = pd.concat([filtered_df, group[either_loop_present]])

    if do_print:
        print("For all the other ones, kept at least %.2f of all the sequences" % fraction_filtered_for_printing)

    return filtered_df


def apply_4_consecutive_filters(
    inp_df,
    max_probability_disrupted = 0.2,
    min_probablity_rescued = 0.7,
    do_allow_coefficient = True,
    disrupted_coefficient = 1,
    rescued_coefficient = 1,
    negligible_probability = 0.3,
    minimal_mutation_size = 2,
    fraction_filtered_for_printing = 0.000001,
    do_print = False):
    crit_1 = filter_non_disruptive_mutations(
        inp_df,
        max_probability_disrupted,
        do_allow_coefficient,
        disrupted_coefficient,
        fraction_filtered_for_printing,
        do_print)

    crit_2 = filter_bad_rescue_mutations(
        crit_1,
        min_probablity_rescued,
        do_allow_coefficient,
        rescued_coefficient,
        fraction_filtered_for_printing,
        do_print)

    crit_3 = filter_double_disruptive_mutations(
        crit_2,
        negligible_probability,
        fraction_filtered_for_printing,
        do_print)

    crit_4 = filter_mutations_by_length(
        crit_3,
        minimal_mutation_size)
    return crit_4


def apply_filters_sort_by_nt_content_take_top_1(
    inp_df,
    max_probability_disrupted = 0.2,
    min_probablity_rescued = 0.7,
    do_allow_coefficient = True,
    disrupted_coefficient = 1,
    rescued_coefficient = 1,
    negligible_probability = 0.3,
    minimal_mutation_size = 2,
    fraction_filtered_for_printing = 0.000001,
    n_top_mutations = 1):
    filtered_df = apply_4_consecutive_filters(
                inp_df,
                max_probability_disrupted,
                min_probablity_rescued,
                do_allow_coefficient,
                disrupted_coefficient,
                rescued_coefficient,
                negligible_probability,
                minimal_mutation_size,
                fraction_filtered_for_printing)
    out_df = rank_mutations_by_change_of_GC_content(filtered_df,
                                                                               n_top_mutations)
    return out_df

def merge_with_loop_selection_df(mut_filtered_df,
                                 loop_selection_df):
    out_df = pd.merge(loop_selection_df, mut_filtered_df,
                      how='inner',
                      left_on=['fr_name', 'source', 'loop'],
                      right_on=['fr_name', 'source', 'changed_loop'])
    out_df.drop(['gene_name_x'], axis=1, inplace=True)
    out_df = out_df.rename({'gene_name_y': 'gene_name'}, axis=1)
    out_df = out_df.sort_values(by=['gene_name', 'fr_name'])
    return out_df


def filter_mutations_by_length(mut_prob_df,
                               minimal_mutation_size = 2
                               ):
    def calculate_mutation_size(x):
        counter = 0
        for i, k in zip(x['original'],x['mutated']):
            if i != k:
                counter += 1
        return counter

    filtered_df = mut_prob_df.copy()
    filtered_df['mutation_size'] = filtered_df.apply(lambda x:
                                                     calculate_mutation_size(x),
                                                     axis = 1)
    filtered_df = filtered_df[filtered_df['mutation_size'] >= minimal_mutation_size]

    return filtered_df


def get_dict_of_groups(inp_df,
                       grouping_variables):
    indices_dict = inp_df.groupby(["gene_name", "changed_loop", "mut_side"]).groups
    dfs_dict = {x : inp_df.loc[indices_dict[x]] for x in indices_dict}
    return dfs_dict


def process_one_mut_selection(curr_mut_df,
                              collected_df,
                              mutations_necessary,
                              seed = None,
                              do_print = False
                              ):
    if not seed is None:
        np.random.seed(seed)

    # this list of mutations is actually empty: skip
    if curr_mut_df.shape[0] == 0:
        if do_print:
            print("Empty dataframe")
    # this list of mutations has enough: choose enough randomly
    elif (mutations_necessary - collected_df.shape[0]) < curr_mut_df.shape[0]:
        chosen_index = np.random.choice(curr_mut_df.index.tolist(),
                                        size = (mutations_necessary - collected_df.shape[0]),
                                        replace = False)
        collected_df = collected_df.append(curr_mut_df.loc[chosen_index],
                                           ignore_index = True)
        if do_print:
            print("Randomly chosen %d mutations" % len(chosen_index))
    # this list of mutations does not have enough: take all
    else:
        collected_df = collected_df.append(curr_mut_df, ignore_index = True)
        if do_print:
            print("Took all the available %d mutations" % curr_mut_df.shape[0])

    return collected_df


def reformat_mutation_selection_dataframe(selection_df):
    reformatted_df = pd.DataFrame(columns = ["gene_name", "fr_name",
                                             "source", "sequence",
                                             "seq_name"])
    for index, row in selection_df.iterrows():
        base_of_row = {"gene_name" : row["gene_name"],
                       "fr_name" : row["fr_name"],
                       "source" : row["source"]}
        mutated_row =  copy.deepcopy(base_of_row)
        rescued_row = copy.deepcopy(base_of_row)
        mutated_row["sequence"] = row['mutated']
        rescued_row["sequence"] = row['rescued']

        if row['mut_name'] == "original":
            mutated_row["seq_name"] = "%s_original" % (row['gene_name'])
            reformatted_df = reformatted_df.append(mutated_row, ignore_index = True)
        else:
            mutated_row["seq_name"] = "%s_%s_mutated" % (row['gene_name'], row['mut_name'])
            rescued_row["seq_name"] = "%s_%s_rescued" % (row['gene_name'], row['mut_name'])
            reformatted_df = reformatted_df.append(mutated_row, ignore_index = True)
            reformatted_df = reformatted_df.append(rescued_row, ignore_index = True)
    return reformatted_df


def check_for_restriction_sites(inp_df,
                                seq_column_name,
                                sites_dict = {"MluI" : "ACGCGT", "PacI" : "TTAATTAA"},
                                do_print = False):
    are_there_rest_sites = pd.Series(index = inp_df.index).fillna(True)
    if do_print:
        print("Number of sequences with restriction sites in them: %d" %
              (are_there_rest_sites.shape[0] - are_there_rest_sites.sum()))
    for enzyme in sites_dict:
        curr_series = inp_df.apply(lambda x:
                                   not (sites_dict[enzyme] in x[seq_column_name]),
                                   axis = 1)
        are_there_rest_sites = are_there_rest_sites & curr_series
        if do_print:
            print("Number of sequences with restriction sites in them: %d" %
                (are_there_rest_sites.shape[0] - are_there_rest_sites.sum()))
    return are_there_rest_sites


def filter_mutations_with_defined_parameters(inp_df,
                                   group_id,
                                  parameters_dict):
    curr_par_dict = parameters_dict[group_id]
    inp_df_crit_1 = filter_non_disruptive_mutations(
                                    inp_df,
                                    max_probability_disrupted = curr_par_dict["max_probability_disrupted"],
                                    coefficient = curr_par_dict["coefficient"])

    inp_df_crit_2 = filter_bad_rescue_mutations(
                                    inp_df_crit_1,
                                    min_probablity_rescued = curr_par_dict["min_probablity_rescued"],
                                    coefficient = curr_par_dict["coefficient"])

    inp_df_crit_3 = filter_double_disruptive_mutations(
                                    inp_df_crit_2,
                                    negligible_probability = curr_par_dict["negligible_probability"])
    return inp_df_crit_3


def calculate_mutations_nt_content(inp_df,
                                   original_nt_content):
    out_df = inp_df.copy()
    out_df['mutated_nt_content'] = out_df.apply(lambda x:
                                                sw_utils.calculate_nt_content(x['mutated']),
                                                axis=1)
    out_df['rescued_nt_content'] = out_df.apply(lambda x:
                                                sw_utils.calculate_nt_content(x['rescued']),
                                                axis=1)

    out_df['mutated_nt_content_difference'] = out_df.apply(lambda x:
                                                           np.sum(
                                                               np.abs(x["mutated_nt_content"] - original_nt_content)),
                                                           axis=1)
    out_df['rescued_nt_content_difference'] = out_df.apply(lambda x:
                                                           np.sum(
                                                               np.abs(x["rescued_nt_content"] - original_nt_content)),
                                                           axis=1)

    out_df['both_nt_content_difference'] = out_df['mutated_nt_content_difference'] + \
                                           out_df['rescued_nt_content_difference']

    return out_df


def check_proportion_top_bottom_10_perc(samples_organization_dict,
                                        values_dict,
                                        fraction = 0.1):
    conf_1_sample = samples_organization_dict["conformation 1"]
    conf_2_sample = samples_organization_dict["conformation 2"]
    conf_1_values = values_dict[conf_1_sample]
    conf_2_values = values_dict[conf_2_sample]
    df_sample_1 = pd.DataFrame({"label" : [1] * conf_1_values.shape[0],
                                "value" : conf_1_values})
    df_sample_2 = pd.DataFrame({"label" : [2] * conf_2_values.shape[0],
                                "value" : conf_2_values})
    combined_df = pd.concat([df_sample_1, df_sample_2])
    combined_df = combined_df.sort_values(by = "value")
    n_in_fraction = round(combined_df.shape[0] * fraction)
    background_fraction_conf_1 = conf_1_values.shape[0] / combined_df.shape[0]
    bottom_fraction = combined_df.iloc[ : n_in_fraction]
    top_fraction = combined_df.iloc[(combined_df.shape[0] - n_in_fraction) : ]
    bottom_fraction_conf_1 = bottom_fraction[bottom_fraction['label'] == 1].shape[0] / bottom_fraction.shape[0]
    top_fraction_conf_1 = top_fraction[top_fraction['label'] == 1].shape[0] / top_fraction.shape[0]

    string_to_print = "default fraction of conformation 1: %.2f\n" % (background_fraction_conf_1)
    string_to_print += "fraction of conformation 1 in the bottom 10 percent: %.2f\n" % (bottom_fraction_conf_1)
    string_to_print += "fraction of conformation 1 in the top 10 percent: %.2f\n" % (top_fraction_conf_1)
    print(string_to_print)


def subset_center_mass_five_perturbations(source_row,
                                     df_with_center_mass_values,
                                     pert_types_order
                                     ):
    general_sub_df = df_with_center_mass_values[(df_with_center_mass_values["fr_name"] == source_row["fr_name"]) &
                                    # or original
                                    ((df_with_center_mass_values["source"] == source_row["source"]) |
                                     (df_with_center_mass_values["source"] == "orig"))].copy()
    general_sub_df.index = general_sub_df['pert_type']
    # to make sure that if a perturbation is missing, an error doesn't get raised
    # for example, in case of 3 perurbations only
    present_mutations_order = [x for x in pert_types_order if x in general_sub_df.index]
    general_sub_df_ordered = general_sub_df.loc[present_mutations_order]
    mass_center_dict = pd.Series(general_sub_df_ordered['dna_mass_center'].values, index=general_sub_df_ordered['pert_type']).to_dict()
    return mass_center_dict


def calculate_shifts_for_perturbations(mutations_df,
                                               counts_df,
                                               perturbations_order,
                                               averaging_method = 'median'
                                               ):
    shifts_df = pd.DataFrame(index=mutations_df.index,
                           columns=['fr_name', 'source', 'gene_name', 'average_shift'])

    shifts_df["fr_name"] = mutations_df['fr_name']
    shifts_df["source"] = mutations_df['source']
    shifts_df["gene_name"] = mutations_df['gene_name']

    for index, row in mutations_df.iterrows():
        center_of_mass_dict = subset_center_mass_five_perturbations(
                                              row,
                                              counts_df,
                                              perturbations_order
                                              )
        orig_center_mass = center_of_mass_dict['orig']
        the_other_centers = np.array([center_of_mass_dict[x] for x in center_of_mass_dict if x != 'orig'])
        shifts_compared_to_orig = the_other_centers - orig_center_mass
        shift_abs_values = np.abs(shifts_compared_to_orig)
        if averaging_method == 'median':
            avg_shift_value = np.median(shift_abs_values)
        elif averaging_method == 'mean':
            avg_shift_value = np.mean(shift_abs_values)
        else:
            print("unknown averaging method")
            sys.error(1)
        shifts_df.at[index, 'average_shift'] = avg_shift_value
    return shifts_df


def calculate_correlation_for_single_perturbations(counts_df,
                                                   replicate_1_column_names,
                                                   replicate_2_column_names
                                                   ):
    correlations_list = []
    for index, row in counts_df.iterrows():
        orig_row = counts_df[(counts_df["fr_name"] == row["fr_name"]) &
                                  (counts_df["source"] == "orig")]
        current_1 = row[replicate_1_column_names].to_numpy().flatten()
        current_2 = row[replicate_2_column_names].to_numpy().flatten()
        orig_1 = orig_row[replicate_1_column_names].to_numpy().flatten()
        orig_2 = orig_row[replicate_2_column_names].to_numpy().flatten()
        r_1, pval_1 = pearsonr(current_1, orig_1)
        r_2, pval_2 = pearsonr(current_2, orig_2)
        average_correlation = (r_1 + r_2) / 2
        correlations_list.append(average_correlation)
    out_df = counts_df.copy()
    out_df['shift'] = correlations_list
    return out_df


def subset_counts_five_perturbations_v2(source_row,
                                     df_with_values,
                                     pert_types_order,
                                     replicate_1_column_names,
                                     replicate_2_column_names,
                                     return_subset_order = False
                                     ):
    general_sub_df = df_with_values[(df_with_values["fr_name"] == source_row["fr_name"]) &
                                    # or original
                                    ((df_with_values["source"] == source_row["source"]) |
                                     (df_with_values["source"] == "orig"))].copy()
    general_sub_df.index = general_sub_df['pert_type']
    # to make sure that if a perturbation is missing, an error doesn't get raised
    # for example, in case of 3 perurbations only
    present_mutations_order = [x for x in pert_types_order if x in general_sub_df.index]
    general_sub_df_ordered = general_sub_df.loc[present_mutations_order]
    rep_1_counts = general_sub_df_ordered[replicate_1_column_names]
    rep_2_counts = general_sub_df_ordered[replicate_2_column_names]
    return rep_1_counts, rep_2_counts, present_mutations_order


def rank_mutations_by_change_of_GC_content(mut_prob_df,
                                           n_top_mutations = 1):
    grouped = mut_prob_df.groupby(["fr_name", "source", "changed_loop"])
    
    def return_top_mutation_from_group(group):
        curr_df = group.copy()
        curr_df['original_nt_content'] = curr_df.apply(lambda x:
                                              sw_utils.calculate_nt_content(x['original']),
                                              axis = 1)
        curr_df['mutated_nt_content'] = curr_df.apply(lambda x:
                                              sw_utils.calculate_nt_content(x['mutated']),
                                              axis = 1)
        curr_df['rescued_nt_content'] = curr_df.apply(lambda x:
                                              sw_utils.calculate_nt_content(x['rescued']),
                                              axis = 1)
        curr_df['mutated_nt_content_difference'] = curr_df.apply(lambda x:
                                                               np.sum(
                                                                   np.abs(
                                                                       x["mutated_nt_content"] - x['original_nt_content'])),
                                                               axis=1)
        curr_df['rescued_nt_content_difference'] = curr_df.apply(lambda x:
                                                               np.sum(
                                                                   np.abs(
                                                                       x["rescued_nt_content"] - x['original_nt_content'])),
                                                               axis=1)

        curr_df['both_nt_content_difference'] = curr_df['mutated_nt_content_difference'] + \
                                               curr_df['rescued_nt_content_difference']
        curr_df = curr_df.sort_values(by = 'both_nt_content_difference')
        top_current_mutation = curr_df.head(n = n_top_mutations)
        return top_current_mutation

    aggregated_top_mutations = grouped.aggregate(return_top_mutation_from_group)
    # clear out the multilevel index
    aggregated_top_mutations = aggregated_top_mutations.droplevel([0, 1, 2], axis=0).reset_index().drop(['index'], axis=1)

    return aggregated_top_mutations


def get_difference_two_mut_selections(new_selection,
                                      old_selection):
    identifiers_df = pd.concat([new_selection[['fr_name', 'source', 'loop']],
               old_selection[['fr_name', 'source', 'loop']]]).drop_duplicates(keep=False)
    out_df = pd.DataFrame(columns = new_selection.columns)

    for index, row in identifiers_df.iterrows():
        curr_full_row = new_selection[
        (new_selection['fr_name'] == row['fr_name']) &
        (new_selection['source'] == row['source']) &
        (new_selection['loop'] == row['loop'])
        ]
        if curr_full_row.shape[0] == 0:
            continue
        assert curr_full_row.shape[0] == 1
        out_df = pd.concat([out_df, curr_full_row], axis = 0)
    return out_df

               