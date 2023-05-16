import numpy as np
import pandas as pd
import sys

sys.path.append('/switchfinder/')
from SwitchFinder.classes import loop, base_pair_array
import SwitchFinder.perturbations as perturbations
import SwitchFinder.dotbracket_comparisons as dotbracket_comparisons
import SwitchFinder.mibp_comparisons as mibp_comparisons


def are_there_no_changing_intervals(inp_list_of_fragments, inp_collection):
    no_ch_int_dict = {}
    for fr_name in inp_list_of_fragments:
        if fr_name not in inp_collection.body_dict:
            no_ch_int_dict[fr_name] = False
            continue
        current_changing_intervals = inp_collection.body_dict[fr_name].get_differential_intervals()
        assert current_changing_intervals.ndim == 2
        curr_no_ch_intervals = current_changing_intervals.shape[0] == 0
        no_ch_int_dict[fr_name] = curr_no_ch_intervals
    no_ch_int_series = pd.Series(no_ch_int_dict)
    return no_ch_int_series


def are_conformations_too_similar(inp_list_of_fragments,
                                inp_collection,
                                min_distance = 8,
                                do_print = False
                                ):
    # min distance of 8 allows switches like this:
    # ........((((......))))....
    # ..((((............))))....
    conform_distances_dict = {}
    for fr_name in inp_list_of_fragments:
        if not fr_name in inp_collection.body_dict:
            conform_distances_dict[fr_name] = 100
            continue
        major_second_distance = dotbracket_comparisons.difference_between_two_conformations(
                                            inp_collection.body_dict[fr_name].major_conf,
                                            inp_collection.body_dict[fr_name].second_conf
                                            )
        conform_distances_dict[fr_name] = major_second_distance
        if do_print:
            if major_second_distance < min_distance:
                print(fr_name)
                print(inp_collection.body_dict[fr_name].major_conf.get_changing_string())
                print(inp_collection.body_dict[fr_name].second_conf.get_changing_string())
                print()

    conform_distances_series = pd.Series(conform_distances_dict)
    are_too_similar_series = conform_distances_series < min_distance
    return are_too_similar_series


def do_conformations_correspond_to_constraint_loops(inp_list_of_fragments,
                                                    mibp_confl_dict,
                                                    inp_collection,
                                                    min_fraction_bps_matching = 0.80,
                                                    do_print = False
                                                    ):
    does_match_dict = {}
    for fr_name in inp_list_of_fragments:
        # check if it's present in collection and in foling
        if not fr_name in inp_collection.body_dict or not fr_name in mibp_confl_dict:
            does_match_dict[fr_name] = True
            continue

        # initialize loops from MIBP output
        first_loop = loop()
        first_loop.initialize_from_dict(mibp_confl_dict[fr_name]['major'])
        second_loop = loop()
        second_loop.initialize_from_dict(mibp_confl_dict[fr_name]['second'])

        # initialize bp objects from the folding output
        major_folding_bp = base_pair_array(inp_collection.body_dict[fr_name].major_conf.string)
        second_folding_bp = base_pair_array(inp_collection.body_dict[fr_name].second_conf.string)

        # do the comparison
        major_does_folding_match, _ = perturbations.does_loop_match_folding(
                                        first_loop,
                                        major_folding_bp,
                                        min_fraction_bps_matching=min_fraction_bps_matching)
        second_does_folding_match, _ = perturbations.does_loop_match_folding(
                                        second_loop,
                                        second_folding_bp,
                                        min_fraction_bps_matching=min_fraction_bps_matching)

        do_both_foldings_match = major_does_folding_match & second_does_folding_match
        does_match_dict[fr_name] = do_both_foldings_match

        if do_print:
            print(fr_name)
            first_loop.print()
            second_loop.print()
            print("major expected")
            print(first_loop.get_dot_bracket_string(np.zeros(major_folding_bp.array.shape[0])))
            print("second expected")
            print(second_loop.get_dot_bracket_string(np.zeros(second_folding_bp.array.shape[0])))
            print("major folded")
            print(major_folding_bp.array_to_string())
            print("second folded")
            print(second_folding_bp.array_to_string())

    does_match_series = pd.Series(does_match_dict)
    return does_match_series


def is_switch_in_the_end(inp_list_of_fragments,
                        mibp_confl_dict,
                        sequences_dict,
                        distance_to_chop_off,
                        do_print = False):
    is_in_the_end_dict = {}
    for fr_name in inp_list_of_fragments:
        # check if it's present in collection and in foling
        if not fr_name in mibp_confl_dict:
            is_in_the_end_dict[fr_name] = False
            continue

        # initialize loops from MIBP output
        first_loop = loop()
        first_loop.initialize_from_dict(mibp_confl_dict[fr_name]['major'])
        second_loop = loop()
        second_loop.initialize_from_dict(mibp_confl_dict[fr_name]['second'])

        min_value = min(first_loop.left_end, second_loop.left_end)
        max_value = max(first_loop.right_end, second_loop.right_end)
        available_on_the_left = min_value
        available_on_the_right = len(sequences_dict[fr_name]) - max_value
        available_total = available_on_the_left + available_on_the_right

        if available_total < distance_to_chop_off:
            is_in_the_end_dict[fr_name] = True
            if do_print:
                print(fr_name)
                print(first_loop.get_dot_bracket_string(np.zeros(len(sequences_dict[fr_name]))))
                print(second_loop.get_dot_bracket_string(np.zeros(len(sequences_dict[fr_name]))))
        else:
            is_in_the_end_dict[fr_name] = False
    is_in_the_end_series = pd.Series(is_in_the_end_dict)
    return is_in_the_end_series


def do_alt_conformations_match(inp_list_of_fragments,
                               collection_of_interest,
                               collection_naive):
    alt_conformation_match_dict = {}

    for fr_name in inp_list_of_fragments:
        if (not fr_name in collection_of_interest.body_dict) or \
            (not fr_name in collection_naive.body_dict):
            alt_conformation_match_dict[fr_name] = False
            continue
        alt_conformation_match_dict[fr_name] = mibp_comparisons.are_fragments_the_same(
                                    collection_of_interest.body_dict[fr_name],
                                    collection_naive.body_dict[fr_name])
    alt_conformation_match_series = pd.Series(alt_conformation_match_dict)
    return alt_conformation_match_series


def compile_denying_criteria_into_one_column(denying_df):
    pass_series = denying_df['in_mibp'] & \
        ~denying_df['too_much_overlap'] & \
        ~denying_df['no_perturb'] & \
        ~denying_df['no_changing_intervals'] & \
        ~denying_df['conformations_too_similar'] & \
        denying_df['correct_folding'] & \
        ~denying_df['switch_on_the_ends']
    return pass_series


def count_number_of_perturbations(perturb_df):
    pert_number_series = (perturb_df['strengthen_passed'].astype(int) + perturb_df['weaken_passed'].astype(int)) * 2
    return pert_number_series


def remove_duplicates(features_df, inp_selection_df,
                      doubled_indication_col,
                      first_selection_col,
                      second_selection_col,
                      first_major_col,
                      second_major_col,
                      first_reciprocal_col,
                      second_reciprocal_col,
                      ):
    selection_df = inp_selection_df.copy()
    doubled_ones = features_df[doubled_indication_col] & inp_selection_df[first_selection_col] & inp_selection_df[second_selection_col]
    processed = pd.Series(index = doubled_ones.index)
    processed = processed.fillna(False)
    # consecutive rules
    # if both are major, remove the second one
    # if one is major, remove the other one
    # if both are reciprocal, remove the second one
    # if one of the two is reciprocal, remove the other one
    # if none are reciprocal, remove the second one
    both_major = doubled_ones & features_df[first_major_col] & features_df[second_major_col]
    selection_df.loc[both_major, second_selection_col] = False
    processed = processed | both_major

    first_major = doubled_ones & ~processed & features_df[first_major_col] & ~features_df[second_major_col]
    selection_df.loc[first_major, second_selection_col] = False
    processed = processed | first_major

    second_major = doubled_ones & ~processed & ~features_df[first_major_col] & features_df[second_major_col]
    selection_df.loc[second_major, first_selection_col] = False
    processed = processed | second_major

    both_reciprocal = doubled_ones & ~processed & features_df[first_reciprocal_col] & features_df[second_reciprocal_col]
    selection_df.loc[both_reciprocal, second_selection_col] = False
    processed = processed | both_reciprocal

    first_reciprocal = doubled_ones & ~processed & features_df[first_reciprocal_col] & ~features_df[second_reciprocal_col]
    selection_df.loc[first_reciprocal, second_selection_col] = False
    processed = processed | first_reciprocal

    second_reciprocal = doubled_ones & ~processed & ~features_df[first_reciprocal_col] & features_df[second_reciprocal_col]
    selection_df.loc[second_reciprocal, first_selection_col] = False
    processed = processed | second_reciprocal

    selection_df.loc[doubled_ones & ~processed, second_selection_col] = False

    assert (selection_df.loc[doubled_ones, [first_selection_col, second_selection_col]].sum(axis=1) == 1).all
    return selection_df


def number_of_perturbations_in_selection(features_df, inp_selection_df):
    n_pert_rep_1 = inp_selection_df['rep_1'] * features_df['n_perturb_rep_1']
    n_pert_rep_2 = inp_selection_df['rep_2'] * features_df['n_perturb_rep_2']
    n_pert_avg = inp_selection_df['avg'] * features_df['n_perturb_avg']
    fragments_chosen = inp_selection_df.sum(axis=1) > 0
    n_pert_total_series = n_pert_rep_1 + n_pert_rep_2 + n_pert_avg + fragments_chosen
    n_pert_total = n_pert_total_series.sum()
    return n_pert_total


def remove_similar_switches(features_df, inp_selection_df):
    # priorities: I trust either single replicate more than the averaged profile
    # however, if one of them is reciprocal I prefer it
    selection_df = inp_selection_df.copy()

    selection_df = remove_duplicates(features_df, selection_df,
                                      doubled_indication_col = 'r1_avg_same_prediction',
                                      first_selection_col = 'rep_1',
                                      second_selection_col = 'avg',
                                      first_major_col = 'r1_major',
                                      second_major_col = 'avg_major',
                                      first_reciprocal_col = 'r1_reciprocal',
                                      second_reciprocal_col = 'avg_reciprocal',
                                      )

    selection_df = remove_duplicates(features_df, selection_df,
                                     doubled_indication_col='r2_avg_same_prediction',
                                     first_selection_col='rep_2',
                                     second_selection_col='avg',
                                     first_major_col='r2_major',
                                     second_major_col='avg_major',
                                     first_reciprocal_col='r2_reciprocal',
                                     second_reciprocal_col='avg_reciprocal',
                                     )

    selection_df = remove_duplicates(features_df, selection_df,
                                     doubled_indication_col='reps_same_prediction',
                                     first_selection_col='rep_1',
                                     second_selection_col='rep_2',
                                     first_major_col='r1_major',
                                     second_major_col='r2_major',
                                     first_reciprocal_col='r1_reciprocal',
                                     second_reciprocal_col='r2_reciprocal',
                                     )

    selection_df = remove_duplicates(features_df, selection_df,
                                      doubled_indication_col = 'small_distance_rep_1_avg',
                                      first_selection_col = 'rep_1',
                                      second_selection_col = 'avg',
                                      first_major_col = 'r1_major',
                                      second_major_col = 'avg_major',
                                      first_reciprocal_col = 'r1_reciprocal',
                                      second_reciprocal_col = 'avg_reciprocal',
                                      )

    selection_df = remove_duplicates(features_df, selection_df,
                                     doubled_indication_col='small_distance_rep_2_avg',
                                     first_selection_col='rep_2',
                                     second_selection_col='avg',
                                     first_major_col='r2_major',
                                     second_major_col='avg_major',
                                     first_reciprocal_col='r2_reciprocal',
                                     second_reciprocal_col='avg_reciprocal',
                                     )

    selection_df = remove_duplicates(features_df, selection_df,
                                     doubled_indication_col='small_distance_rep_1_rep_2',
                                     first_selection_col='rep_1',
                                     second_selection_col='rep_2',
                                     first_major_col='r1_major',
                                     second_major_col='r2_major',
                                     first_reciprocal_col='r1_reciprocal',
                                     second_reciprocal_col='r2_reciprocal',
                                     )

    return selection_df


def initiate_acceptable_masks(features_df):
    # create selection dataframe
    selection_df = pd.DataFrame(index = features_df.index,
                                columns = ['rep_1','rep_2','avg'])
    selection_df = selection_df.fillna(False)

    # which fragments are acceptable by coverage and denying criteria
    rep_1_acceptable = features_df['covered_1'] & features_df['passes_DMS_1']
    rep_2_acceptable = features_df['covered_2'] & features_df['passes_DMS_2']
    average_acceptable = (features_df['covered_1'] | features_df['covered_2']) & features_df['passes_avg']
    return selection_df, rep_1_acceptable, rep_2_acceptable, average_acceptable


def classify_fragments_by_switch_categories_main_group(features_df):
    selection_df, rep_1_acceptable, rep_2_acceptable, average_acceptable = initiate_acceptable_masks(features_df)
    # which fragments pass repiprocal test
    rep_1_recipr = features_df['r1_reciprocal'] & rep_1_acceptable
    rep_2_recipr = features_df['r2_reciprocal'] & rep_2_acceptable
    avg_recipr = features_df['avg_reciprocal'] & average_acceptable

    # all fragments that pass reciprocal test and are acceptable should be selected
    selection_df.loc[rep_1_recipr, 'rep_1'] = True
    selection_df.loc[rep_2_recipr, 'rep_2'] = True
    selection_df.loc[avg_recipr, 'avg'] = True

    # for fragments where only one replicate is reciprocal and average is acceptable: add the average
    selection_df.loc[(rep_1_recipr & ~rep_2_recipr & average_acceptable), 'avg'] = True
    selection_df.loc[(~rep_1_recipr & rep_2_recipr & average_acceptable), 'avg'] = True

    selection_df = remove_similar_switches(features_df, selection_df)
    n_perturbations = number_of_perturbations_in_selection(features_df, selection_df)

    return selection_df, n_perturbations


def classify_fragments_by_switch_categories_other_groups(features_df,
                                                         mask_of_interest,
                                                         mask_already_processed):
    selection_df, rep_1_acceptable, rep_2_acceptable, average_acceptable = initiate_acceptable_masks(features_df)
    rep_1_acceptable = rep_1_acceptable & ~mask_already_processed & mask_of_interest
    rep_2_acceptable = rep_2_acceptable & ~mask_already_processed & mask_of_interest
    average_acceptable = average_acceptable & ~mask_already_processed & mask_of_interest

    selection_df.loc[average_acceptable, 'avg'] = True
    selection_df.loc[(rep_1_acceptable & ~average_acceptable & ~rep_2_acceptable), 'rep_1'] = True
    selection_df.loc[(~rep_1_acceptable & ~average_acceptable & rep_2_acceptable), 'rep_2'] = True
    selection_df.loc[(rep_1_acceptable & ~average_acceptable & rep_2_acceptable), 'rep_2'] = True
    n_perturbations = number_of_perturbations_in_selection(features_df, selection_df)

    return selection_df, n_perturbations


def classify_fragments_by_switch_categories_other_groups_one_rep(features_df,
                                                         mask_of_interest,
                                                         mask_already_processed,
                                                         rep):
    selection_df, rep_1_acceptable, rep_2_acceptable, average_acceptable = initiate_acceptable_masks(features_df)

    if rep == 'rep_1':
        rep_1_acceptable = rep_1_acceptable & ~mask_already_processed & mask_of_interest
        selection_df.loc[rep_1_acceptable, 'rep_1'] = True
    elif rep == 'rep_2':
        rep_2_acceptable = rep_2_acceptable & ~mask_already_processed & mask_of_interest
        selection_df.loc[rep_2_acceptable, 'rep_2'] = True
    else:
        print("wrong replicate!")
        sys.exit(1)

    n_perturbations = number_of_perturbations_in_selection(features_df, selection_df)

    return selection_df, n_perturbations


def combine_two_selection_dataframes(df_1, df_2):
    assert (df_1.index == df_2.index).all()
    columns_list = ['rep_1', 'rep_2', 'avg']
    new_selection_df = pd.DataFrame(index=df_1.index,
                                columns=columns_list)
    for col in columns_list:
        new_selection_df[col] = df_1[col] | df_2[col]

    return new_selection_df


def combine_selection_dataframes(features_df,
                                 list_of_selection_dfs):
    selection_df = pd.DataFrame(index=features_df.index,
                                columns=['rep_1', 'rep_2', 'avg'])
    selection_df = selection_df.fillna(False)

    for curr_df in list_of_selection_dfs:
        selection_df = combine_two_selection_dataframes(selection_df, curr_df)
    return selection_df


def choose_dots_within_std_ellipse(inp_df,
                                   columns_of_interest,
                                   number_of_stds):
    assert len(columns_of_interest) == 2
    assert len(number_of_stds) == 2

    distance_metric = np.zeros(inp_df.shape[0])
    for i, col_name in enumerate(columns_of_interest):
        curr_column_median = inp_df[col_name].median()
        curr_column_std = inp_df[col_name].std()
        curr_number_of_stds = number_of_stds[i]
        curr_radius = curr_number_of_stds * curr_column_std
        curr_column_np = inp_df[col_name].to_numpy()
        current_additive = ((curr_column_np - curr_column_median) ** 2) / (curr_radius ** 2)
        distance_metric += current_additive

    return distance_metric <= 1


def get_perturbation_sequences_from_selection(selection_df,
                                              pert_dfs_dict,
                                              orig_seq_dict):
    out_df = pd.DataFrame(columns = ['fr_name', 'source', 'pert_type', 'sequence'])
    pert_strength_types_list = ["major_strengthen_pert",
                                "second_strengthen_pert"]
    pert_weak_types_list = ["second_weaken_pert",
                                "major_weaken_pert"]

    for fr_name in selection_df.index:
        if selection_df.loc[fr_name].sum() == 0:
            continue
        orig_sequence = orig_seq_dict[fr_name]
        orig_element = {'fr_name' : fr_name,
                         'source': 'orig',
                         'pert_type' : 'orig',
                         'sequence' : orig_sequence}
        out_df = out_df.append(orig_element, ignore_index=True)

        for col in selection_df:
            if not selection_df.loc[fr_name, col]:
                continue
            curr_perturb_row = pert_dfs_dict[col].loc[fr_name]
            if curr_perturb_row['strengthen_passed']:
                for strengthen_type in pert_strength_types_list:
                    current_sequence = curr_perturb_row[strengthen_type]
                    curr_perturbation = {'fr_name' : fr_name,
                                         'source': col,
                                         'pert_type' : strengthen_type,
                                         'sequence' : current_sequence}
                    out_df = out_df.append(curr_perturbation, ignore_index=True)

            if curr_perturb_row['weaken_passed']:
                for strengthen_type in pert_weak_types_list:
                    current_sequence = curr_perturb_row[strengthen_type]
                    curr_perturbation = {'fr_name' : fr_name,
                                         'source': col,
                                         'pert_type' : strengthen_type,
                                         'sequence' : current_sequence}
                    out_df = out_df.append(curr_perturbation, ignore_index=True)

    return out_df


def get_leftmost_and_rightmost_changing_coords(
                    pert_df_row,
                    conflicts_general_dict,
                    values_to_fill_in = (100, 0),
                    do_print = False):
    fr_name = pert_df_row['fr_name']
    source = pert_df_row['source']

    # if it's original fragment, we don't have one conflict item for it casue it depends on which sources have been chosen
    if source not in conflicts_general_dict:
        return values_to_fill_in[0], values_to_fill_in[1]

    conflict_dict = conflicts_general_dict[source][fr_name]

    # initialize loops
    first_loop = loop()
    first_loop.initialize_from_dict(conflict_dict['major'])
    second_loop = loop()
    second_loop.initialize_from_dict(conflict_dict['second'])

    min_value = min(first_loop.left_end, second_loop.left_end)
    max_value = max(first_loop.right_end, second_loop.right_end)

    return min_value, max_value


def get_shortening_coordinates_for_a_fragment(changing_coords_row,
                                              expected_length = 186,
                                              to_chop_off = 12):
    left_coord = 0
    right_coord = changing_coords_row['length']

    # make sure fragment isn't already short
    if changing_coords_row['length'] > expected_length - to_chop_off:
        # which side has more space available
        start_with_left = changing_coords_row['available_left'] > changing_coords_row['available_right']
        if start_with_left:
            if changing_coords_row['available_left'] >= to_chop_off:
                left_coord = to_chop_off
            else:
                left_coord = changing_coords_row['min_changing']
                right_coord = changing_coords_row['length'] - (to_chop_off - changing_coords_row['min_changing'])
        else:
            if changing_coords_row['available_right'] >= to_chop_off:
                right_coord = changing_coords_row['length'] - to_chop_off
            else:
                right_coord = changing_coords_row['max_changing']
                left_coord = to_chop_off - changing_coords_row['available_right']

    left_coord = int(left_coord)
    right_coord = int(right_coord)

    assert (right_coord - left_coord) <= (expected_length - to_chop_off)
    assert changing_coords_row['min_changing'] >= left_coord
    assert changing_coords_row['max_changing'] <= right_coord

    return left_coord, right_coord


def cut_a_sequence(sequence,
                  fr_name,
                  cut_df):
    cut_sub_df = cut_df[cut_df['fr_name'] == fr_name]
    assert cut_sub_df.shape[0] == 1
    left_coord = cut_sub_df['left_coord'].item()
    right_coord = cut_sub_df['right_coord'].item()
    subsequence = sequence[left_coord : right_coord]

    return subsequence
