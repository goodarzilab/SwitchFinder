import numpy as np
import os
import sys

import glob_vars
import folding_api

def tuple_to_int(inp_tuple):
    temp_list = [int(x) for x in inp_tuple]
    return tuple(temp_list)


# from here https://stackoverflow.com/questions/43600878/merging-overlapping-intervals
# note that the array has to be sorted!
def merge_intervals(intervals):
    starts = intervals[:,0]
    ends = np.maximum.accumulate(intervals[:,1])
    valid = np.zeros(len(intervals) + 1, dtype=bool)
    valid[0] = True
    valid[-1] = True
    valid[1:-1] = starts[1:] >= ends[:-1]
    return np.vstack((starts[:][valid[:-1]], ends[:][valid[1:]])).T


def get_confl_loops_dict(inp_dict):
    out_dict = {}
    for fr in inp_dict:
        out_dict[fr] = folding_api.get_two_major_conflicting_loops(inp_dict[fr])
    return out_dict

def array_to_string(inp_array):
    string_array = [glob_vars._char_to_nt_mapping[x] for x in inp_array]
    return "".join(string_array)


def string_to_array(inp_string):
    out_array = np.zeros(len(inp_string), dtype=np.uint8)
    for i in range(len(inp_string)):
        out_array[i] = glob_vars._nt_to_char_mapping[inp_string[i]]

    return out_array


def convert_conflict_dict_to_intervals(conflict_dict):
    intervals_array = np.zeros((4, 2))
    intervals_array[0, 0] = conflict_dict['major']['forward'][0]
    intervals_array[0, 1] = conflict_dict['major']['forward'][1]
    intervals_array[1, 0] = conflict_dict['major']['backward'][0]
    intervals_array[1, 1] = conflict_dict['major']['backward'][1]
    intervals_array[2, 0] = conflict_dict['second']['forward'][0]
    intervals_array[2, 1] = conflict_dict['second']['forward'][1]
    intervals_array[3, 0] = conflict_dict['second']['backward'][0]
    intervals_array[3, 1] = conflict_dict['second']['backward'][1]

    # make it integer
    intervals_array = intervals_array.astype(np.int)
    # make them open-ended
    intervals_array[:, 1] += 1
    return intervals_array


def check_if_switch_in_the_end(switches_df,
                              conflicts_dict,
                              sequences_dict,
                              distance_to_chop_off):
    out_df = switches_df.copy()
    counter = 0
    for index, row in switches_df.iterrows():
        for col in switches_df.columns:
            source_of_mibp = row[col]
            if source_of_mibp == '':
                continue
            current_conflict = conflicts_dict[source_of_mibp][row.name]
            assert len(current_conflict) > 0

            intervals_array = convert_conflict_dict_to_intervals(current_conflict)
            intervals_array.sort(axis=0)
            merged_intervals_array = merge_intervals(intervals_array)
            min_value = np.min(merged_intervals_array)
            max_value = np.max(merged_intervals_array)
            available_on_the_left = min_value
            available_on_the_right = len(sequences_dict[row.name]) - max_value
            available_total = available_on_the_left + available_on_the_right

            if available_total < distance_to_chop_off:
                out_df.at[row.name, col] = ''
                counter += 1
    print("Total switches filtered out: %d" % counter)
    return out_df


def get_GC_content(inp_seq):
    return (inp_seq.count("C") + inp_seq.count("G")) / len(inp_seq)


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def geom_mean_1d(inp_array):
    x = np.log(inp_array[inp_array > 0]).sum() / inp_array.shape[0]
    return np.exp(x)


def calculate_nt_content(seq_string):
    counts = np.zeros(4, dtype = np.float)
    for i, nt in enumerate(glob_vars.nt_list):
        nt_string = glob_vars._char_to_nt_mapping[nt]
        counts[i] = seq_string.count(nt_string)

    counts = counts / counts.sum()

    return counts


# from here https://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-whilst-preserving-order
def remove_duplicates_list_preserve_order(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]