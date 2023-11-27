import numpy as np
import textwrap

import sys


import SwitchFinder.glob_vars as glob_vars


def get_stretch_coordinates(current_stretch):
    first_stretch_beginning = current_stretch[0][0]
    first_stretch_end = current_stretch[-1][0]
    second_stretch_beginning = current_stretch[-1][1]
    second_stretch_end = current_stretch[0][1]
    return [first_stretch_beginning, first_stretch_end, \
           second_stretch_beginning, second_stretch_end]


def list_all_stretches_numpy(current_structure):
    positions_stack = []
    all_stretches_list = []
    current_stretch = []
    for i in range(0, len(current_structure)):
        letter = current_structure[i]
        if letter == '.':
            if len(positions_stack) > 0:
                if positions_stack[-1] != '.':
                    positions_stack.append('.')

        elif letter == '(':
            positions_stack.append(i)

        elif letter == ')':
            if len(positions_stack) == 0:
                print("Error! A closing bracket without opening bracket found!\n\n\n")
                break

            if positions_stack[-1] == '.':
                positions_stack.pop()
                if len(current_stretch) != 0:
                    current_stretch = current_stretch[::-1]
                    current_stretch_array = get_stretch_coordinates(current_stretch)
                    all_stretches_list.append(current_stretch_array)
                current_stretch = []

            first_base = positions_stack[-1]
            second_base = i

            basepair = (first_base, second_base)
            current_stretch.append(basepair)
            positions_stack.pop()

        else:
            print("Error! An unknown character was found!\n\n\n")
            break

    if len(current_stretch) != 0:
        current_stretch = current_stretch[::-1]
        current_stretch_array = get_stretch_coordinates(current_stretch)
        all_stretches_list.append(current_stretch_array)

    all_stretches_list.sort(key=lambda x: x[1])
    all_stretches_np = np.array(all_stretches_list, dtype=int)

    return all_stretches_np


def make_sure_there_are_no_other_stretches(current_structure, current_stretch, next_stretch):
    beg_left_interval = current_stretch[1] + 1
    end_left_interval = next_stretch[0]
    if beg_left_interval > end_left_interval:
        return False
    left_interval_string = current_structure[beg_left_interval: end_left_interval]

    beg_right_interval = next_stretch[3] + 1
    end_right_interval = current_stretch[2]
    right_interval_string = current_structure[beg_right_interval: end_right_interval]

    if beg_right_interval > end_right_interval:
        return False

    if '(' not in left_interval_string and \
                    ')' not in left_interval_string and \
                    '(' not in right_interval_string and \
                    ')' not in right_interval_string:
        return True


def combine_two_stretches(stretches_np_loc, i):
    loc_list = []
    for k in range(i):
        loc_list.append(stretches_np_loc[k])
    cur_str = stretches_np_loc[i]
    next_str = stretches_np_loc[i + 1]
    combined_stretch = [cur_str[0], next_str[1], next_str[2], cur_str[3]]
    loc_list.append(combined_stretch)
    for k in range(i + 2, stretches_np_loc.shape[0]):
        loc_list.append(stretches_np_loc[k])

    return np.array(loc_list, dtype=int)


def find_parts_of_the_same_stretch(current_structure, stretches_np_inp, MAX_SPACING):
    stretches_np_loc = stretches_np_inp.copy()
    #     print(stretches_np_loc)
    #     for k in stretches_np_loc:
    #         print(current_structure[k[0] : k[3] + 1])

    number_of_stretches = stretches_np_loc.shape[0]

    any_common_stretches_found = False

    for i in range(number_of_stretches - 1):
        current_stretch = stretches_np_loc[i]
        next_stretch = stretches_np_loc[i + 1]
        distance_between_left = next_stretch[0] - current_stretch[1] - 1
        distance_between_right = current_stretch[2] - next_stretch[3] - 1
        if (distance_between_left <= MAX_SPACING) and \
                (distance_between_right <= MAX_SPACING):
            if make_sure_there_are_no_other_stretches(current_structure, current_stretch, next_stretch):
                #                 print('\n------\n')
                #                 print(current_structure[current_stretch[0] : current_stretch[3] + 1])
                #                 print(current_structure[next_stretch[0] : next_stretch[3] + 1])
                #                 print(i)
                #                 print('\n------\n')


                any_common_stretches_found = True

                stretches_np_loc = combine_two_stretches(stretches_np_loc, i)
                break


    if any_common_stretches_found:
        stretches_np_loc = find_parts_of_the_same_stretch(current_structure, stretches_np_loc, MAX_SPACING)
    else:
        return stretches_np_loc

    return stretches_np_loc


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def remove_changing_parts_both_conformations(spaced_stretches_1_inp, spaced_stretches_2_inp,
                                                inp_string_1, inp_string_2,
                                             MAXIMAL_OVERHANG):
    indices_same_1 = []
    indices_same_2 = []

    #     print(spaced_stretches_1_inp)
    #     print(spaced_stretches_2_inp)
    #     print()

    counter_1 = 0

    for i in spaced_stretches_1_inp:
        counter_1 += 1
        counter_2 = 0

        for k in spaced_stretches_2_inp:
            counter_2 += 1

            overlap_left = getOverlap((i[0], i[1] + 1), (k[0], k[1] + 1))
            overlap_right = getOverlap((i[2], i[3] + 1), (k[2], k[3] + 1))

            if (overlap_left == 0) or (overlap_right == 0):
                continue

            left_stem1_length = i[1] + 1 - i[0]
            left_stem2_length = k[1] + 1 - k[0]
            max_left_stem_length = max(left_stem1_length, left_stem2_length)

            right_stem1_length = i[3] + 1 - i[2]
            right_stem2_length = k[3] + 1 - k[2]
            max_right_stem_length = max(right_stem1_length, right_stem2_length)

            if (max_left_stem_length - overlap_left > MAXIMAL_OVERHANG) or \
                    (max_right_stem_length - overlap_right > MAXIMAL_OVERHANG):
                # print(i)
                # print(k)
                # print(inp_string_1[i[0] : i[3] + 1])
                # print(inp_string_2[k[0] : k[3] + 1])
                # print(textwrap.fill(inp_string_1, width=91))
                # print(textwrap.fill(inp_string_2, width=91))
                # print(overlap_left, overlap_right)
                # print()
                continue

            indices_same_1.append(counter_1 - 1)
            indices_same_2.append(counter_2 - 1)


    indices_same_1 = sorted(list(set(indices_same_1)))
    indices_same_2 = sorted(list(set(indices_same_2)))

    all_indices_1 = list(range(spaced_stretches_1_inp.shape[0]))
    all_indices_2 = list(range(spaced_stretches_2_inp.shape[0]))

    indices_different_1 = sorted(list(set(all_indices_1).difference(set(indices_same_1))))
    indices_different_2 = sorted(list(set(all_indices_2).difference(set(indices_same_2))))

    unchanged_stretches_1 = spaced_stretches_1_inp[indices_same_1]
    changed_stretches_1 = spaced_stretches_1_inp[indices_different_1]
    unchanged_stretches_2 = spaced_stretches_2_inp[indices_same_2]
    changed_stretches_2 = spaced_stretches_2_inp[indices_different_2]

    return unchanged_stretches_1, changed_stretches_1, \
           unchanged_stretches_2, changed_stretches_2


def string_and_stretches_to_print(inp_string, stretches_loc):
    structure_list_original = list(inp_string)
    structure_list_copy = ['.'] * len(structure_list_original)

    for i in stretches_loc:
        structure_list_copy[i[0]: i[1] + 1] = structure_list_original[i[0]: i[1] + 1]
        structure_list_copy[i[2]: i[3] + 1] = structure_list_original[i[2]: i[3] + 1]

    return ''.join(structure_list_copy)


def remove_short_changing_stems(stretches_loc, MIN_STEM_LENGHT):
    indices_short = []

    counter = 0

    for i in stretches_loc:
        counter += 1
        left_stem_length = i[1] + 1 - i[0]
        right_stem_length = i[3] + 1 - i[2]

        min_stem_length = min(left_stem_length, right_stem_length)

        if min_stem_length <= MIN_STEM_LENGHT:
            indices_short.append(counter - 1)

    all_indices = list(range(stretches_loc.shape[0]))
    indices_not_short = sorted(list(set(all_indices).difference(set(indices_short))))
    filtered_stretches_loc = stretches_loc[indices_not_short]

    return filtered_stretches_loc


def check_if_the_stem_is_outmost(stem_of_interest,
                                 other_stems_1,
                                 other_stems_2,
                                 other_stems_3):
    if other_stems_1.shape[0] == 0:
        return False

    is_outmost = True
    for np_stems in [other_stems_1, other_stems_2, other_stems_3]:
        if np_stems.shape[0] == 0:
            continue
        if (stem_of_interest[1] >= np_stems[0][0]) or \
                (stem_of_interest[2] <= np_stems[-1][-1]):
            is_outmost = False

    return is_outmost


def remove_one_outer_stem(unchanged_stretches,
                          changed_stretches_1,
                          changed_stretches_2):
    # print(unchanged_stretches)
    the_outmost_stem = unchanged_stretches[0]

    if not check_if_the_stem_is_outmost(the_outmost_stem,
                                        unchanged_stretches[1:],
                                        changed_stretches_1,
                                        changed_stretches_2):
        return unchanged_stretches
    return unchanged_stretches[1:]


def remove_all_outer_stems(unchanged_stretches,
                           changed_stretches_1,
                           changed_stretches_2):
    current_stretches = unchanged_stretches.copy()
    trimmed_stretches = remove_one_outer_stem(current_stretches,
                                              changed_stretches_1,
                                              changed_stretches_2)

    while trimmed_stretches.shape[0] != current_stretches.shape[0]:
        current_stretches = trimmed_stretches.copy()
        trimmed_stretches = remove_one_outer_stem(current_stretches,
                                                  changed_stretches_1,
                                                  changed_stretches_2)
    return trimmed_stretches


def two_stretches_difference(first_array, second_array):
    indices_to_remove = []

    for i in range(first_array.shape[0]):
        for k in range(second_array.shape[0]):
            if (first_array[i] == second_array[k]).all():
                indices_to_remove.append(i)

    indices_to_remove = sorted(list(set(indices_to_remove)))

    all_indices = list(range(first_array.shape[0]))
    indices_not_short = sorted(list(set(all_indices).difference(set(indices_to_remove))))
    filtered_stretches_loc = first_array[indices_not_short]

    return filtered_stretches_loc


def convert_string_to_pairing_states(inp_string):
    pairing_states_array = np.zeros(len(inp_string), dtype=np.uint8)
    for i, ch in enumerate(inp_string):
        pairing_states_array[i] = glob_vars._char_to_structure[ch]
    return pairing_states_array


# def convert_string_to_pairing_states_with_stem_ends(
#                                 pairing_states_array,
#                                 all_stretches):
#     extended_pairing_states_array = pairing_states_array.copy()
#     for stem_tuple in all_stretches:
#         left_stem_beginning = stem_tuple[0]
#         left_stem_end = stem_tuple[1]
#         right_stem_beginning = stem_tuple[2]
#         right_stem_end = stem_tuple[3]
#
#         for coord in [left_stem_beginning, left_stem_end,
#                       right_stem_beginning, right_stem_end]:
#             extended_pairing_states_array[coord] = glob_vars._end_of_stem
#
#     return extended_pairing_states_array



def convert_string_to_numpy(inp_string):
    numpy_array = np.zeros(len(inp_string), dtype=np.uint8)
    for i, ch in enumerate(inp_string):
        numpy_array[i] = glob_vars._char_to_numpy[ch]
    return numpy_array


def identify_differentially_paired_positions(pairing_states_1, pairing_states_2):
    unpaired_in_first = pairing_states_1 == glob_vars._loop
    unpaired_in_second = pairing_states_2 == glob_vars._loop
    all_differential_positions = pairing_states_1 != pairing_states_2
    diff_unpaired_1 = np.logical_and(all_differential_positions, unpaired_in_first)
    diff_unpaired_2 = np.logical_and(all_differential_positions, unpaired_in_second)

    return (diff_unpaired_1, diff_unpaired_2)


def identify_constantly_paired_positions(pairing_states_1, pairing_states_2):
    all_constant_positions = pairing_states_1 == pairing_states_2
    constantly_paired_positions = np.logical_and(all_constant_positions, pairing_states_1 == glob_vars._stem)
    constantly_unpaired_positions = np.logical_and(all_constant_positions, pairing_states_1 == glob_vars._loop)
    return constantly_paired_positions, constantly_unpaired_positions


def difference_between_two_conformations(conf_1, conf_2):
    if conf_1.numpy is None:
        conf_1.convert_to_numpy()
    if conf_2.numpy is None:
        conf_2.convert_to_numpy()
    assert conf_1.numpy.shape[0] == conf_2.numpy.shape[0]
    difference = (conf_1.numpy != conf_2.numpy).sum()
    return difference
