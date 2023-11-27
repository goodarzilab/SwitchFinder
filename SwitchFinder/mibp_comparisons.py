import os, sys
import numpy as np
import pandas as pd
import copy

current_script_path = sys.argv[0]
package_home_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..'))
if package_home_path not in sys.path:
    sys.path.append(package_home_path)

from sw_finder.wrappers import fold_alternative_consensus_structures_parralel as fold_script
from . import glob_vars, utils, dotbracket_comparisons


def compare_two_intervals(tuple_1, tuple_2, allowed_distance = 1):
    assert len(tuple_1) == 2
    assert len(tuple_2) == 2
    difference_beginning = abs(tuple_1[0] - tuple_2[0])
    difference_end = abs(tuple_1[1] - tuple_2[1])
    are_two_intervals_equal = False
    if difference_beginning <= allowed_distance and difference_end <= allowed_distance:
        are_two_intervals_equal = True
    return are_two_intervals_equal


def compare_two_loops(loop_1, loop_2,
                      take_MI_into_account = False,
                      allowed_distance = 1):
    are_loops_equal = True
    if not compare_two_intervals(loop_1['forward'], loop_2['forward'],
                                 allowed_distance = allowed_distance):
        are_loops_equal = False
    if not compare_two_intervals(loop_1['backward'], loop_2['backward'],
                                 allowed_distance = allowed_distance):
        are_loops_equal = False
    if take_MI_into_account:
        if not compare_two_intervals(loop_1['average_MI'], loop_2['average_MI'],
                                     allowed_distance=allowed_distance):
            are_loops_equal = False
    return are_loops_equal


def compare_if_major_loops_are_the_same(entry_1, entry_2,
                                        take_MI_into_account = False,
                                        allowed_distance = 1):
    entry_1_major_loop, entry_1_major_loop_MI = fold_script.get_major_loop(entry_1)
    entry_2_major_loop, entry_2_major_loop_MI = fold_script.get_major_loop(entry_2)

    entry_1_second_loop = fold_script.get_max_loop_conflicting_with_the_major_one(
                                                entry_1,
                                                entry_1_major_loop,
                                                entry_1_major_loop_MI)
    entry_2_second_loop = fold_script.get_max_loop_conflicting_with_the_major_one(
                                                entry_2,
                                                entry_2_major_loop,
                                                entry_2_major_loop_MI)

    first_orientation = True
    if not compare_two_loops(entry_1_major_loop, entry_2_major_loop,
                      take_MI_into_account = take_MI_into_account,
                      allowed_distance = allowed_distance):
        first_orientation = False
    if not compare_two_loops(entry_1_second_loop, entry_2_second_loop,
                      take_MI_into_account = take_MI_into_account,
                      allowed_distance=allowed_distance):
        first_orientation = False

    second_orientation = True
    if not compare_two_loops(entry_1_second_loop, entry_2_major_loop,
                      take_MI_into_account = take_MI_into_account,
                      allowed_distance = allowed_distance):
        second_orientation = False
    if not compare_two_loops(entry_1_major_loop, entry_2_second_loop,
                      take_MI_into_account = take_MI_into_account,
                      allowed_distance=allowed_distance):
        second_orientation = False
    are_elements_equal = first_orientation or second_orientation
    return are_elements_equal


def compare_all_the_loops(entry_1, entry_2,
                          take_MI_into_account = False,
                          allowed_distance = 1):
    are_elements_equal = True
    pair_loops_entry_1 = entry_1['pairs_loops']
    pair_loops_entry_2 = entry_2['pairs_loops']
    if len(pair_loops_entry_1) != len(pair_loops_entry_2):
        are_elements_equal = False
        return are_elements_equal

    for num1 in pair_loops_entry_1:
        for num2 in pair_loops_entry_1[num1]:
            current_loop_entry_1 = pair_loops_entry_1[num1][num2]
            current_loop_entry_2 = pair_loops_entry_2[num1][num2]
            if not compare_two_loops(current_loop_entry_1, current_loop_entry_2,
                                     take_MI_into_account=take_MI_into_account,
                                     allowed_distance=allowed_distance):
                are_elements_equal = False
    return are_elements_equal


def compare_if_major_loop_is_in_second_entry(entry_1, entry_2,
                          take_MI_into_account = False,
                          allowed_distance = 1):
    is_major_pair_loop_detected = False
    entry_1_major_loop, entry_1_major_loop_MI = fold_script.get_major_loop(entry_1)
    entry_1_second_loop = fold_script.get_max_loop_conflicting_with_the_major_one(
                                                entry_1,
                                                entry_1_major_loop,
                                                entry_1_major_loop_MI)
    pair_loops_entry_2 = entry_2['pairs_loops']
    for num1 in pair_loops_entry_2:
        for num2 in pair_loops_entry_2[num1]:
            current_loop_entry_1 = pair_loops_entry_2[num1][num2]
            if compare_two_loops(current_loop_entry_1, entry_1_major_loop,
                                     take_MI_into_account=take_MI_into_account,
                                     allowed_distance=allowed_distance):
                # get the complementary loop
                major_loop_numb_1 = num1
                major_loop_numb_2 = num2
                possible_values_num2 = copy.deepcopy(set(pair_loops_entry_2[num1].keys()))
                possible_values_num2.remove(major_loop_numb_2)
                assert len(possible_values_num2) == 1
                possible_second_loop_num2 = list(possible_values_num2)[0]
                possible_second_loop = pair_loops_entry_2[major_loop_numb_1][possible_second_loop_num2]
                if compare_two_loops(possible_second_loop, entry_1_second_loop,
                                     take_MI_into_account=take_MI_into_account,
                                     allowed_distance=allowed_distance):
                    is_major_pair_loop_detected = True
    return is_major_pair_loop_detected




def compare_if_major_loops_are_the_same_dicts(dict_1, dict_2,
                                              take_MI_into_account = False):
    assert len(dict_1) == len(dict_2)
    out_dict = {}
    for el in dict_1:
        out_dict[el] = compare_if_major_loops_are_the_same(dict_1[el], dict_2[el],
                                                           take_MI_into_account = take_MI_into_account)
    return out_dict


def compare_all_the_loops_dicts(dict_1, dict_2,
                                take_MI_into_account = False,
                                allowed_distance = 1):
    assert len(dict_1) == len(dict_2)
    out_dict = {}
    for el in dict_1:
        out_dict[el] = compare_all_the_loops(dict_1[el], dict_2[el],
                                        take_MI_into_account = take_MI_into_account,
                                        allowed_distance=allowed_distance)
    return out_dict


def compare_MIBP_with_and_without_DMS(with_dms_dict,
                                      no_dms_dict,
                                      seq_dict):
    with_without_dms_comparison_dict = {}

    for fr_name in seq_dict:
        curr_dict_element = {}
        curr_dict_element['in_no_dms'] = fr_name in no_dms_dict
        if fr_name in with_dms_dict:
            curr_dict_element['in_dms'] = True
            if curr_dict_element['in_no_dms']:
                curr_dict_element['all_loops'] = compare_all_the_loops(
                                                        no_dms_dict[fr_name],
                                                        with_dms_dict[fr_name],
                                                        take_MI_into_account = False,
                                                        allowed_distance=1)
                curr_dict_element['major_loops'] = compare_if_major_loops_are_the_same(
                                                        no_dms_dict[fr_name],
                                                        with_dms_dict[fr_name],
                                                        take_MI_into_account = False,
                                                        allowed_distance=1)
                curr_dict_element['old_major_loop_in_dms'] = compare_if_major_loop_is_in_second_entry(
                                                        no_dms_dict[fr_name],
                                                        with_dms_dict[fr_name],
                                                        take_MI_into_account = False,
                                                        allowed_distance=1)
                curr_dict_element['dms_major_loop_in_orig'] = compare_if_major_loop_is_in_second_entry(
                                                        with_dms_dict[fr_name],
                                                        no_dms_dict[fr_name],
                                                        take_MI_into_account = False,
                                                        allowed_distance=1)

            else:
                curr_dict_element['all_loops'] = False
                curr_dict_element['major_loops'] = False
                curr_dict_element['old_major_loop_in_dms'] = False
                curr_dict_element['dms_major_loop_in_orig'] = False
        else:
            curr_dict_element['in_dms'] = False
            curr_dict_element['all_loops'] = False
            curr_dict_element['major_loops'] = False
            curr_dict_element['old_major_loop_in_dms'] = False
            curr_dict_element['dms_major_loop_in_orig'] = False
        curr_dict_element['original'] = fr_name.endswith('original')
        with_without_dms_comparison_dict[fr_name] = curr_dict_element

    with_without_dms_comparison_df = pd.DataFrame.from_dict(with_without_dms_comparison_dict,
                                                            orient = 'index')

    return with_without_dms_comparison_df


def apply_reciprocal_criteria(inp_df,
                              col_name_1 = 'dms_major_loop_in_orig',
                              col_name_2 = 'old_major_loop_in_dms',
                              do_return_column = False):
    out_df = inp_df.copy()
    out_column = out_df[col_name_1] & out_df[col_name_2]
    out_df = out_df[out_column]
    if do_return_column:
        return out_column
    else:
        return out_df


def apply_major_criteria(inp_df):
    out_df = inp_df.copy()
    out_df = out_df[out_df['major_loops']]
    return out_df


def print_fraction_of_original(inp_df):
    n_original = inp_df['original'].sum()
    fraction = n_original / inp_df.shape[0]
    print("Fraction of original fragments is %.2f" % (fraction))


def print_number_unique_vs_common_index(df_1, df_2):
    intersection = df_1.index.intersection(df_2.index)
    difference_left = df_1.index.difference(df_2.index)
    difference_right = df_2.index.difference(df_1.index)
    difference_size = len(difference_left) + len(difference_right)
    fraction = len(intersection) / min(df_1.shape[0], df_2.shape[0])
    print("Dataframes' sizes: %d, %d" % (df_1.shape[0], df_2.shape[0]))
    print("Indexes common between two dataframes: %d (fraction %.2f); indexes unique %d" % \
          (len(intersection), fraction, difference_size))


def count_number_of_original_fragments_with_their_counterparts(inp_df):
    original_sub_df = inp_df[inp_df['original']]
    shuffled_sub_df = inp_df[~inp_df['original']]
    original_names = set([x.replace('_original', '') for x in original_sub_df.index])
    shuffled_names = set([x.replace('_shuffling_0', '') for x in shuffled_sub_df.index])
    intersection = original_names.intersection(shuffled_names)
    fraction = len(intersection) / len(original_names)
    print("%d out of %d original fragments (fraction %.2f) among the ones that passed have their counterparts also pass" % \
          (len(intersection), len(original_names), fraction))


def are_structures_the_same(conf_1, conf_2,
                            maximal_distance_between_structures):
    are_the_same = False
    difference = (conf_1.numpy != conf_2.numpy).sum()
    if difference <= maximal_distance_between_structures:
        are_the_same = True
    return are_the_same



def are_fragments_the_same(fr_1, fr_2,
                           maximal_distance_between_structures = 2):
    fr_1.convert_to_numpy()
    fr_2.convert_to_numpy()

    are_conformations_the_same_direct_order = True
    if not are_structures_the_same(fr_1.major_conf, fr_2.major_conf,
                            maximal_distance_between_structures):
        are_conformations_the_same_direct_order = False
    if not are_structures_the_same(fr_1.second_conf, fr_2.second_conf,
                            maximal_distance_between_structures):
        are_conformations_the_same_direct_order = False

    are_conformations_the_same_reverse_order = True
    if not are_structures_the_same(fr_1.major_conf, fr_2.second_conf,
                            maximal_distance_between_structures):
        are_conformations_the_same_reverse_order = False
    if not are_structures_the_same(fr_1.second_conf, fr_2.major_conf,
                            maximal_distance_between_structures):
        are_conformations_the_same_reverse_order = False

    are_fragments_the_same = are_conformations_the_same_direct_order or are_conformations_the_same_reverse_order
    return are_fragments_the_same


def difference_between_folds_from_two_collections(fr,
                                                  collect_1, collect_2,
                                                  negative_value = -1):
    if (not fr in collect_1.body_dict) or (not fr in collect_2.body_dict):
        return negative_value
    major_1_major_2_dist = dotbracket_comparisons.difference_between_two_conformations(
                                    collect_1.body_dict[fr].major_conf,
                                    collect_2.body_dict[fr].major_conf)
    major_1_second_2_dist = dotbracket_comparisons.difference_between_two_conformations(
                                collect_1.body_dict[fr].major_conf,
                                collect_2.body_dict[fr].second_conf)
    major_2_second_1_dist = dotbracket_comparisons.difference_between_two_conformations(
                                    collect_1.body_dict[fr].second_conf,
                                    collect_2.body_dict[fr].major_conf)
    second_1_second_2_dist = dotbracket_comparisons.difference_between_two_conformations(
                                collect_1.body_dict[fr].second_conf,
                                collect_2.body_dict[fr].second_conf)
    direct_dist = (major_1_major_2_dist + second_1_second_2_dist) / 2
    reversed_dist = (major_2_second_1_dist + major_1_second_2_dist) / 2
    return min(direct_dist, reversed_dist)
