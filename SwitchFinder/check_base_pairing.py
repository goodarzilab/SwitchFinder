import numpy as np
import os
import sys

import sys
sys.path.append('/switchfinder/')

import SwitchFinder.glob_vars

def make_reverse_complement(sequence_np):
    out_sequence = np.zeros_like(sequence_np)
    for i in range(sequence_np.shape[0]):
        out_sequence[i] = glob_vars._base_pairing_dict[sequence_np[i]]
    reversed = np.flip(out_sequence)
    return reversed


def is_paired(base1, base2,
              allow_canonical = True):
    is_paired_bool = False
    if (base1 == glob_vars._T and base2 == glob_vars._A) or \
        (base1 == glob_vars._C and base2 == glob_vars._G) or \
        (base1 == glob_vars._G and base2 == glob_vars._C) or \
        (base1 == glob_vars._A and base2 == glob_vars._T):
        is_paired_bool = True
    if allow_canonical:
        if is_paired_bool or \
        (base1 == glob_vars._T and base2 == glob_vars._G) or \
        (base1 == glob_vars._G and base2 == glob_vars._T):
            is_paired_bool = True
    return is_paired_bool


def n_not_paired_bases(sequence_np_1, sequence_np_2,
                       orientation = "reversed",
                       allow_canonical = True):
    sequence_np_2_c = sequence_np_2.copy()
    assert sequence_np_1.shape[0] == sequence_np_2_c.shape[0]
    is_pared_array = np.zeros_like(sequence_np_1, dtype=np.bool)
    if orientation == "reversed":
        sequence_np_2_c = np.flip(sequence_np_2_c)
    for i in range(sequence_np_1.shape[0]):
        is_pared_array[i] = is_paired(sequence_np_1[i], sequence_np_2_c[i],
                                      allow_canonical = allow_canonical)
    number_of_mismatches = np.invert(is_pared_array).sum()
    return number_of_mismatches


# see here https://stackoverflow.com/questions/24885092/finding-the-consecutive-zeros-in-a-numpy-array
def true_runs(inp_array):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    padded = np.zeros(inp_array.shape[0] + 2)
    padded[1:-1] = inp_array
    absdiff = np.abs(np.diff(padded))
    ranges = np.where(absdiff == 1)[0].reshape(-1, 2)
    return ranges


def longest_true_run(inp_array):
    ranges = true_runs(inp_array)
    if ranges.shape[0] == 0:
        return 0
    length = ranges[: , 1] - ranges[: , 0]
    max_length = np.max(length)
    return max_length


def max_base_pairs_in_a_row(sequence_np_1, sequence_np_2,
                       orientation = "reversed",
                       allow_canonical = True):
    sequence_np_2_c = sequence_np_2.copy()
    assert sequence_np_1.shape[0] == sequence_np_2_c.shape[0]
    is_pared_array = np.zeros_like(sequence_np_1, dtype=np.bool)
    if orientation == "reversed":
        sequence_np_2_c = np.flip(sequence_np_2_c)
    for i in range(sequence_np_1.shape[0]):
        is_pared_array[i] = is_paired(sequence_np_1[i], sequence_np_2_c[i],
                                      allow_canonical = allow_canonical)
    longest_run = longest_true_run(is_pared_array)
    return longest_run


def are_segments_paired(sequence_np_1, sequence_np_2,
                        orientation = "reversed",
                        allow_canonical = True):
    number_of_mismatches = n_not_paired_bases(sequence_np_1, sequence_np_2,
                                              orientation = orientation,
                                              allow_canonical = allow_canonical)
    are_paired = number_of_mismatches == 0
    return are_paired


def hamming_distance_two_seqs(sequence_np_1, sequence_np_2):
    return (sequence_np_1 != sequence_np_2).sum()


