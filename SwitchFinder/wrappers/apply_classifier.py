import sys
sys.path.append('/switchfinder/')

import argparse
import SwitchFinder.IO as IO
import numpy as np
import pandas as pd

def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies_filename", help="the fasta file with the target sequences", type=str)
    parser.add_argument("--dataframe_output", help="output dataframe filename", type=str)
    parser.add_argument("--text_output", help="output text filename", type=str)
    parser.add_argument("--text_output_short", help="output text filename", type=str)
    parser.add_argument("--loop_energies_coefficient", help="linear model coefficient for the energies of alternative conformations", type=float)
    parser.add_argument("--barrier_heights_coefficient", help="linear model coefficient for the height of the activation energy barrier", type=float)
    parser.add_argument("--intercept", help="output dataframe filename", type=float)

    parser.set_defaults(
                loop_energies_coefficient = 0.128,
                barrier_heights_coefficient = -0.137,
                intercept = 0.542)
    args = parser.parse_args(raw_args)
    return args


def parse_energies_v7_file(inp_energies_file_loc):
    energies_consensus_dict_loc = {}
    fragments_sequence_dict = {}

    with open(inp_energies_file_loc, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n$$$\n')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            splitted_entry = entry.split('\n')
            fragment_name = splitted_entry[0]
            fragment_sequence = splitted_entry[1]
            fragments_sequence_dict[fragment_name] = fragment_sequence

            energies_string_splitted = splitted_entry[8].split(", ")
            major_loop_energy = float(energies_string_splitted[0].split(': ')[1])
            second_loop_energy = float(energies_string_splitted[1].split(': ')[1])
            common_energy = float(energies_string_splitted[2].split(': ')[1])
            if fragment_name not in energies_consensus_dict_loc:
                energies_consensus_dict_loc[fragment_name] = {}
            energies_consensus_dict_loc[fragment_name]['major_loop'] = major_loop_energy
            energies_consensus_dict_loc[fragment_name]['second_loop'] = second_loop_energy
            energies_consensus_dict_loc[fragment_name]['common'] = common_energy

            energies_consensus_dict_loc[fragment_name]['major_loop_consensus'] = splitted_entry[3]
            energies_consensus_dict_loc[fragment_name]['second_loop_consensus'] = splitted_entry[5]

            if (energies_consensus_dict_loc[fragment_name]['major_loop'] == 0) or \
                    (energies_consensus_dict_loc[fragment_name]['second_loop'] == 0):
                del energies_consensus_dict_loc[fragment_name]

    return energies_consensus_dict_loc, fragments_sequence_dict


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


def get_energies_and_barrier_heights(energies_consensus_dict_loc):
    average_loop_energies_dict = {}
    average_barrier_heights_dict = {}
    for fragment_name in energies_consensus_dict_loc:
        major_loop_energy = energies_consensus_dict_loc[fragment_name]['major_loop']
        second_loop_energy = energies_consensus_dict_loc[fragment_name]['second_loop']
        average_loop_energy = -(major_loop_energy + second_loop_energy) / float(2)
        average_loop_energies_dict[fragment_name] = average_loop_energy

        common_energy = energies_consensus_dict_loc[fragment_name]['common']
        first_barrier = common_energy - major_loop_energy
        second_barrier = common_energy - second_loop_energy
        average_barrier = (first_barrier + second_barrier) / float(2)
        average_barrier_heights_dict[fragment_name] = average_barrier

    return average_loop_energies_dict, average_barrier_heights_dict



def two_dicts_to_dataframe(inp_dict_1, inp_dict_2, inp_name_1, inp_name_2):
    df_1 = pd.DataFrame.from_dict(inp_dict_1, orient='index')
    df_1.columns = [inp_name_1]
    df_2 = pd.DataFrame.from_dict(inp_dict_2, orient='index')
    df_2.columns = [inp_name_2]
    #print(df_1.shape, df_2.shape)
    merged_pd = df_1.merge(df_2, how='inner',
                          left_index=True, right_index=True)
    return merged_pd

def normalize_features(inp_dataframe,
                      column_names = ["loop_energies","barrier_heights"]):
    new_df = inp_dataframe.copy()
    for col in column_names:
        curr_series = inp_dataframe[col]
        mean = curr_series.mean()
        std = curr_series.std()
        new_series = (curr_series - mean) / std
        new_df[col] = new_series
    return new_df

def calculate_score(row,
                loop_energies_coefficient,
                barrier_heights_coefficient,
                intercept_coefficient):
    score = (row["loop_energies"] * loop_energies_coefficient +
             row["barrier_heights"] * barrier_heights_coefficient +
             intercept_coefficient)
    return score


def post_process_dataframe(inp_df):
    inp_df = inp_df.rename({"score" : "Classifier score"}, axis = 1)
    out_df = inp_df[["Classifier score"]]
    out_df = out_df.sort_values(by = "Classifier score", ascending = False)
    return out_df


def classify_fragments(features_pd, args):
    features_pd_norm = normalize_features(features_pd)
    features_pd_norm['score'] = features_pd_norm.apply(lambda x:
                                                       calculate_score(x,
                                                       args.loop_energies_coefficient,
                                                       args.barrier_heights_coefficient,
                                                       args.intercept),
                                                       axis=1)
    scores_df = post_process_dataframe(features_pd_norm)
    return scores_df


def convert_energy_calculations_to_df(inp_energies_v4_filename):
    energies_v4_dict_loc, _ = parse_energies_v7_file(inp_energies_v4_filename)
    loop_energies_dict_loc, barrier_heights_dict_loc = get_energies_and_barrier_heights(energies_v4_dict_loc)

    MI_features_df_loc = two_dicts_to_dataframe(
        loop_energies_dict_loc,
        barrier_heights_dict_loc,
        "loop_energies",
        "barrier_heights")
    return MI_features_df_loc
    #MI_features_df_loc.to_csv(out_df_filename, sep='\t')


def make_sting_to_write(fragment_name, fragment_sequence, free_energies_dict,
                        consensus_structures_dict):

    string_to_write = '%s\n%s\n' % (fragment_name, fragment_sequence)
    string_to_write += 'major_loop\n%s\nsecond_loop\n%s\ncommon\n%s\n' % (consensus_structures_dict['major_loop'], consensus_structures_dict['second_loop'], consensus_structures_dict['common'])
    string_to_write += 'major_loop: %f, second_loop: %f, common: %f\n' % (free_energies_dict['major_loop'], free_energies_dict['second_loop'], free_energies_dict['common'])
    string_to_write += ''
    string_to_write += '$$$\n'
    return string_to_write

def write_fragments_in_order(scores_df, inp_energies_v4_filename, text_output_filename):
    energies_v4_dict_loc, _ = parse_energies_v7_file(inp_energies_v4_filename)
    fragments_consensus_structures_dict, fragments_sequence_dict = read_dotbracket_file(inp_energies_v4_filename)
    with open(text_output_filename, 'w') as wf:
        for name in scores_df.index:
            string_to_write = make_sting_to_write(name,
                                fragments_sequence_dict[name],
                                energies_v4_dict_loc[name],
                                fragments_consensus_structures_dict[name])
            wf.write(string_to_write)

def make_sting_to_write_short(fragment_name, fragment_sequence, free_energies_dict,
                        consensus_structures_dict):

    string_to_write = '%s\n%s\n' % (fragment_name, fragment_sequence)
    string_to_write += 'Conformation 1\n%s\nConformation 2\n%s\nCommon substructures\n%s\n' % (consensus_structures_dict['major_loop'], consensus_structures_dict['second_loop'], consensus_structures_dict['common'])
    string_to_write += '$$$\n'
    return string_to_write

def write_fragments_in_order_short(scores_df, inp_energies_v4_filename, text_output_filename):
    energies_v4_dict_loc, _ = parse_energies_v7_file(inp_energies_v4_filename)
    fragments_consensus_structures_dict, fragments_sequence_dict = read_dotbracket_file(inp_energies_v4_filename)
    with open(text_output_filename, 'w') as wf:
        for name in scores_df.index:
            string_to_write = make_sting_to_write_short(name,
                                fragments_sequence_dict[name],
                                energies_v4_dict_loc[name],
                                fragments_consensus_structures_dict[name])
            wf.write(string_to_write)



def main(raw_args = None):
    args = handler(raw_args)
    features_pd = convert_energy_calculations_to_df(args.energies_filename)
    scores_df = classify_fragments(features_pd, args)
    scores_df.to_csv(args.dataframe_output, sep='\t')
    write_fragments_in_order(scores_df, args.energies_filename, args.text_output)
    write_fragments_in_order_short(scores_df, args.energies_filename, args.text_output_short)



if __name__ == "__main__":
    main()