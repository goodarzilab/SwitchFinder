import os, sys
import numpy as np
import re
import time
import math
from subprocess import PIPE, run, Popen, call, list2cmdline
import multiprocessing
import random
import string

import sys
sys.path.append('/switchfinder/')

import SwitchFinder.glob_vars as glob_vars
import SwitchFinder.utils as utils


def extract_probabilities_matrix_from_ps_file(ps_filename,
                                              seq_length,
                                              return_optimal = False):
    # adopted from here: https://github.com/klavinslab/coral/blob/master/coral/analysis/_structure/viennarna.py
    pattern = 'start of base pair probability data\n(.*)\nshowpage'

    with open(ps_filename, 'r') as rf:
        ps_full_string = rf.read()
    # In the default mode, this matches any character except a newline. If the DOTALL flag has been specified, this matches any character including a newline.
    # see here: https://docs.python.org/2/library/re.html
    res = re.search(pattern, ps_full_string, flags=re.DOTALL)
    if res is not None:
        dotplot_lines = res.group(1).split('\n')
    else:
        print("ps file is empty!")
        sys.exit(1)

    ensemble_probs = np.zeros((seq_length, seq_length))
    optimal_probs = np.zeros((seq_length, seq_length))

    for point in dotplot_lines:
        point_split = point.split(' ')
        # Use zero indexing
        i = int(point_split[0]) - 1
        j = int(point_split[1]) - 1
        sqprob = float(point_split[2])
        probtype = point_split[3]
        if probtype == 'ubox':
            ensemble_probs[i][j] = sqprob ** 2
        else:
            optimal_probs[i][j] = sqprob ** 2

    if return_optimal:
        return ensemble_probs, optimal_probs
    else:
        return ensemble_probs


def extract_probabilities_matrix_RNAstructure(probabilities_filename):
    with open(probabilities_filename, 'r') as rf:
        probs_full_string = rf.read()
    os.remove(probabilities_filename)
    probs_full_array = probs_full_string.split('\n')
    seq_length = int(probs_full_array[0])
    ensemble_probs_log = np.zeros((seq_length, seq_length))

    for line in probs_full_array[2:]:
        if line == '':
            continue
        point_split = line.split('\t')
        # Use zero indexing
        i = int(point_split[0]) - 1
        j = int(point_split[1]) - 1
        log10prob = float(point_split[2])
        ensemble_probs_log[i][j] = log10prob

    ensemble_probs = np.power(10, -1 * ensemble_probs_log)
    ensemble_probs[ensemble_probs_log == 0] = 0

    return ensemble_probs




def write_probabilities_to_txt_file(pfs_filename,
                                    RNAstructure_path
                                              ):
    probability_plot_command = os.path.join(RNAstructure_path, "ProbabilityPlot")
    output_filename = pfs_filename.replace("_partition.txt", "_probability_plot.txt")
    probability_plot_arguments = [probability_plot_command, pfs_filename, output_filename, "--text"]
    result = run(args=probability_plot_arguments, stdout=PIPE, stderr=PIPE)
    # print(result.stderr)
    # print(" ".join(probability_plot_arguments))
    if os.path.isfile(pfs_filename):
        os.remove(pfs_filename)
    return output_filename



# modified from sw_finder.wrappers.fold_alternative_consensus_structures_parralel
# this is for ViennaRNA
def write_all_constraints_files(loops_dict_fr, fragment_name, temp_files_folder,
                                include_random_string=True):
    rand_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
    constraints_filenames = []
    for counter in sorted(list(loops_dict_fr.keys())):
        if include_random_string:
            curr_short_filename = fragment_name + '_loop_%s_constraint_%s_ViennaRNA.txt' % (counter, rand_string)
        else:
            curr_short_filename = fragment_name + '_loop_%s_constraint_ViennaRNA.txt' % (counter)
        current_filename = os.path.join(temp_files_folder, curr_short_filename)
        write_constraint_file_for_loop(loops_dict_fr[counter], current_filename)
        constraints_filenames.append(current_filename)

    return constraints_filenames

def write_constraint_file_for_loop(loop_entry, current_filename):
    forward_beginning = int(loop_entry["forward"][0]) + 1
    forward_end = int(loop_entry["forward"][1]) + 1
    backward_beginning = int(loop_entry["backward"][0]) + 1
    backward_end = int(loop_entry["backward"][1]) + 1

    # print(current_filename)
    with open(current_filename, 'w') as wf:
        loop_length = forward_end - forward_beginning + 1
        string_to_write = "F\t%d\t%d\t%d\tA\n" % (forward_beginning, backward_end, loop_length)
        wf.write(string_to_write)

# this is for RNAstructure
def write_all_constraints_files_RNAstructure(loops_dict_fr, fragment_name, temp_files_folder,
                                             include_random_string=True,
                                             no_constraints = False):
    rand_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
    constraints_filenames = []
    if not no_constraints:
        for counter in sorted(list(loops_dict_fr.keys())):
            if include_random_string:
                curr_short_filename = fragment_name + '_loop_%s_constraint_%s_ViennaRNA.txt' % (counter, rand_string)
            else:
                curr_short_filename = fragment_name + '_loop_%s_constraint_ViennaRNA.txt' % (counter)
            current_filename = os.path.join(temp_files_folder, curr_short_filename)
            write_constraint_file_for_loop_RNAstructure(loops_dict_fr[counter], current_filename)
            constraints_filenames.append(current_filename)
    else:
        counter = 0
        if include_random_string:
            curr_short_filename = fragment_name + '_loop_%s_constraint_%s_ViennaRNA.txt' % (counter, rand_string)
        else:
            curr_short_filename = fragment_name + '_loop_%s_constraint_ViennaRNA.txt' % (counter)
        current_filename = os.path.join(temp_files_folder, curr_short_filename)
        constraints_filenames.append(current_filename)

    return constraints_filenames


# this is for RNAstructure
# based on this format description https://rna.urmc.rochester.edu/Text/File_Formats.html#Constraint
def write_constraint_file_for_loop_RNAstructure(loop_entry, current_filename):
    constraints_prefix = "DS:\n-1\nSS:\n-1\nMod:\n-1\nPairs:\n"
    constraints_suffix = "\n-1 -1\nFMN:\n-1\nForbids:\n-1 -1\n"

    forward_beginning = int(loop_entry["forward"][0]) + 1
    forward_end = int(loop_entry["forward"][1]) + 1
    backward_beginning = int(loop_entry["backward"][0]) + 1
    backward_end = int(loop_entry["backward"][1]) + 1
    loop_length = forward_end - forward_beginning + 1

    strings_array = []
    for i in range(loop_length):
        forward_nt = forward_beginning + i
        backward_nt = backward_end - i
        strings_array.append("%d %d" % (forward_nt, backward_nt))
    pairs_string = "\n".join(strings_array)

    string_to_write = "%s%s%s" % (constraints_prefix, pairs_string, constraints_suffix)
    with open(current_filename, 'w') as wf:
        wf.write(string_to_write)


# modified from sw_finder.wrappers.fold_alternative_consensus_structures_parralel
def write_fragment_sequence_file(fragment_name, dict_entry, temp_files_folder):
    rand_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
    current_filename = os.path.join(temp_files_folder, "%s_%s_sequence.fa" % (fragment_name, rand_string))
    with open(current_filename, 'w') as wf:
        wf.write(">%s\n%s\n" % (fragment_name, dict_entry))
    return current_filename


def write_shape_constraints_ViennaRNA(fragment_name, shape_array, temp_files_folder):
    pass
    # shape_filename = os.path.join(temp_files_folder, fragment_name + "_shape.dat")
    # strings_to_write_list = []
    # for i in range(shape_array.shape[0]):
    #     k = i + 1 # move to 1-based coordinate system
    #     strings_to_write_list.append("%d\t%.5f" % (k, shape_array[i]))
    #
    # string_to_write = "\n".join(strings_to_write_list)
    # with open(shape_filename, 'w') as wf:
    #     wf.write(string_to_write)
    #
    # return shape_filename


def write_shape_constraints_RNAstructure(fragment_name, shape_array, temp_files_folder,
                                         include_random_string = True):
    rand_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
    shape_filename = os.path.join(temp_files_folder, "%s_%s_shape.dat" % (fragment_name, rand_string))
    if not include_random_string:
        shape_filename = os.path.join(temp_files_folder, "%s_shape.dat" % (fragment_name))
    strings_to_write_list = []
    for i in range(shape_array.shape[0]):
        k = i + 1 # move to 1-based coordinate system
        strings_to_write_list.append("%d\t%.5f" % (k, shape_array[i]))

    string_to_write = "\n".join(strings_to_write_list)
    with open(shape_filename, 'w') as wf:
        wf.write(string_to_write)

    return shape_filename


def launch_folding_get_probabilities_ViennaRNA(fragment_name,
                                     sequence,
                                     loops_dict_fr,
                                     shape_array,
                                     temp_files_folder,
                                     no_constraints = False,
                                     RNAfold_path = '/rumi/shams/khorms/programs/ViennaRNA-2.4.14/bin/bin/RNAfold',
                                     ):
    constraints_filenames = write_all_constraints_files(loops_dict_fr, fragment_name, temp_files_folder)
    sequence_filename = write_fragment_sequence_file(fragment_name, sequence, temp_files_folder)
    if not shape_array is None:
        shape_filename = write_shape_constraints_ViennaRNA(fragment_name, shape_array, temp_files_folder)

    matrices_list = []
    for constr_filename in constraints_filenames:
        if RNAfold_path is None:
            RNAfold_path = 'RNAfold'
        RNAfold_arguments = [RNAfold_path, '-p', '--bppmThreshold', '0']
        if not no_constraints:
            RNAfold_arguments.append('--constraint=%s' % constr_filename)
        if not shape_array is None:
            RNAfold_arguments.append('--shape')
            RNAfold_arguments.append(shape_filename)
        RNAfold_arguments.append(sequence_filename)
        result = run(args=RNAfold_arguments, stdout=PIPE, stderr=PIPE, cwd=temp_files_folder)
        dot_ps_file = os.path.join(temp_files_folder, '%s_dp.ps' % fragment_name)
        ensemble_probs = extract_probabilities_matrix_from_ps_file(dot_ps_file, len(sequence))
        matrices_list.append(ensemble_probs)
    if no_constraints:
        return matrices_list[0]
    return tuple(matrices_list)


def launch_folding_get_probabilities_RNAstructure(fragment_name,
                                     sequence,
                                     loops_dict_fr,
                                     shape_array,
                                     temp_files_folder,
                                     no_constraints = False,
                                     RNAstructure_path = '/rumi/shams/khorms/programs/RNAstructure',
                                     ):
    RNAstructure_parameters_data_path = os.path.join(RNAstructure_path, 'data_tables')
    RNAstructure_path = os.path.join(RNAstructure_path, 'exe')
    # RNAstructure won't work without DATAPATH environmental variable set to the right folder
    # I can't set it with subprocess.Popen for reasons described here: https://stackoverflow.com/questions/3929319/subprocess-module-errors-with-export-in-python-on-linux
    os.putenv("DATAPATH", RNAstructure_parameters_data_path)
    constraints_filenames = write_all_constraints_files_RNAstructure(loops_dict_fr, fragment_name,
                                                                     temp_files_folder,
                                                                     include_random_string = True,
                                                                     no_constraints = no_constraints)
    sequence_filename = write_fragment_sequence_file(fragment_name, sequence, temp_files_folder)
    if not shape_array is None:
        shape_filename = write_shape_constraints_RNAstructure(fragment_name, shape_array, temp_files_folder)

    matrices_list = []
    for constr_filename in constraints_filenames:
        pfs_filename = constr_filename.replace(".txt", "_partition.txt")
        partition_command = os.path.join(RNAstructure_path, "partition")
        partition_arguments = [partition_command, sequence_filename, pfs_filename]
        if not no_constraints:
            partition_arguments += ['--constraint', constr_filename]

        if not shape_array is None:
            partition_arguments.append('--SHAPE')
            partition_arguments.append(shape_filename)

        result = run(args=partition_arguments, stdout=PIPE, stderr=PIPE, cwd=temp_files_folder)
        # print(" ".join(partition_arguments))
        # print(result.stderr)
        probabilities_filename = write_probabilities_to_txt_file(pfs_filename, RNAstructure_path)
        ensemble_probs = extract_probabilities_matrix_RNAstructure(probabilities_filename)
        matrices_list.append(ensemble_probs)
        if os.path.isfile(constr_filename):
            os.remove(constr_filename)
    if not shape_array is None:
        os.remove(shape_filename)
    os.remove(sequence_filename)

    if no_constraints:
        return matrices_list[0]
    return tuple(matrices_list)


def difference_normalized_by_second_matrix(matrix_1, matrix_2):
    difference = np.sum(np.abs(matrix_1 - matrix_2))
    normalize_by = matrix_2.sum()
    return difference / normalize_by


def get_major_loop(init_pars_dict_entry):
    major_loop = {}
    major_loop_MI = 0
    for num1 in init_pars_dict_entry['pairs_loops']:
        for num2 in init_pars_dict_entry['pairs_loops'][num1]:
            current_loop = init_pars_dict_entry['pairs_loops'][num1][num2]
            current_loop_MI = current_loop['average_MI']
            if current_loop_MI > major_loop_MI:
                major_loop_MI = current_loop_MI
                major_loop = current_loop
    return major_loop, major_loop_MI


def get_characteristic_string(current_loop):
    characteristic_string = "%s\t%s\t%f" % ("\t".join([str(x) for x in current_loop['forward']]),
                            "\t".join([str(x) for x in current_loop['backward']]),
                            current_loop['average_MI'])
    return characteristic_string


def get_max_loop_conflicting_with_the_major_one(init_pars_dict_entry, major_loop, major_loop_MI):
    second_loop = {}
    second_loop_MI = 0
    major_loop_characteristic_string = get_characteristic_string(major_loop)

    for num1 in init_pars_dict_entry['pairs_loops']:
        characteristic_strings_dict = {}
        for num2 in init_pars_dict_entry['pairs_loops'][num1]:
            current_loop = init_pars_dict_entry['pairs_loops'][num1][num2]
            characteristic_string = get_characteristic_string(current_loop)
            characteristic_strings_dict[characteristic_string] = current_loop
        if major_loop_characteristic_string in characteristic_strings_dict:
            del characteristic_strings_dict[major_loop_characteristic_string]
            keys_list = sorted(list(characteristic_strings_dict.keys()))
            current_loop = characteristic_strings_dict[keys_list[0]]
            current_loop_MI = current_loop['average_MI']
            if len(keys_list) != 1:
                print("Error!")
            if current_loop_MI > second_loop_MI:
                second_loop = current_loop
                second_loop_MI = current_loop_MI

    return second_loop


def get_all_base_probabilities_matrices_one_fragment(inp_arguments):
    fr_name = inp_arguments[0]
    seq = inp_arguments[1]
    MIBP_entry = inp_arguments[2]
    shape_profile_1 = inp_arguments[3]
    shape_profile_2 = inp_arguments[4]
    temp_files_folder = inp_arguments[5]

    curr_major_loop, curr_major_loop_MI = get_major_loop(MIBP_entry)
    curr_second_loop = get_max_loop_conflicting_with_the_major_one(
        MIBP_entry, curr_major_loop, curr_major_loop_MI)
    curr_loops_dict = {'major': curr_major_loop, 'second': curr_second_loop}
    out_dict = {}

    m_no_constr = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                seq,
                                                                curr_loops_dict,
                                                                shape_array=None,
                                                                temp_files_folder=temp_files_folder,
                                                                no_constraints=True
                                                                )

    m_shape_1_no_constr = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                        seq,
                                                                        curr_loops_dict,
                                                                        shape_array=shape_profile_1,
                                                                        temp_files_folder=temp_files_folder,
                                                                        no_constraints=True
                                                                        )

    m_shape_2_no_constr = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                        seq,
                                                                        curr_loops_dict,
                                                                        shape_array=shape_profile_2,
                                                                        temp_files_folder=temp_files_folder,
                                                                        no_constraints=True
                                                                        )

    m_no_shape_tuple = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                     seq,
                                                                     curr_loops_dict,
                                                                     shape_array=None,
                                                                     temp_files_folder=temp_files_folder
                                                                     )

    m_shape_1_tuple = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                    seq,
                                                                    curr_loops_dict,
                                                                    shape_array=shape_profile_1,
                                                                    temp_files_folder=temp_files_folder
                                                                    )

    m_shape_2_tuple = launch_folding_get_probabilities_RNAstructure(fr_name,
                                                                    seq,
                                                                    curr_loops_dict,
                                                                    shape_array=shape_profile_2,
                                                                    temp_files_folder=temp_files_folder
                                                                    )
    out_dict['name'] = fr_name
    out_dict['unconstr'] = m_no_constr
    out_dict['DMS1_unconstr'] = m_shape_1_no_constr
    out_dict['DMS2_unconstr'] = m_shape_2_no_constr
    out_dict['major_loop'] = m_no_shape_tuple[0]
    out_dict['second_loop'] = m_no_shape_tuple[1]
    out_dict['DMS1_major_loop'] = m_shape_1_tuple[0]
    out_dict['DMS1_second_loop'] = m_shape_1_tuple[1]
    out_dict['DMS2_major_loop'] = m_shape_2_tuple[0]
    out_dict['DMS2_second_loop'] = m_shape_2_tuple[1]

    return out_dict



def prepare_inputs_for_parallel_base_prob_calculation(
                               seq_dict, MIBP_dict,
                               shape_1_dict, shape_2_dict,
                               temp_files_folder):
    all_fragment_names = sorted(list(seq_dict.keys()))
    arguments_list = []

    for c, fr_name in enumerate(all_fragment_names):
        if fr_name not in MIBP_dict:
            continue
        current_arguments = []
        current_arguments.append(fr_name)
        current_arguments.append(seq_dict[fr_name])
        current_arguments.append(MIBP_dict[fr_name])
        current_arguments.append(shape_1_dict[fr_name])
        current_arguments.append(shape_2_dict[fr_name])
        current_arguments.append(temp_files_folder)
        current_arguments_tuple = tuple(current_arguments)
        arguments_list.append(current_arguments_tuple)

    return arguments_list


def calculate_base_pair_probabilities_many_fragments(
                               seq_dict, MIBP_dict,
                               shape_1_dict, shape_2_dict,
                               temp_files_folder,
                               do_print = True,
                               n_processes = 10,
                               how_often_print = 10
                                ):
    arguments_list = prepare_inputs_for_parallel_base_prob_calculation(
                               seq_dict, MIBP_dict,
                               shape_1_dict, shape_2_dict,
                               temp_files_folder)
    pool = multiprocessing.Pool(n_processes)
    out_dict = {}
    i = 0
    tic = time.time()
    for dict_element in pool.imap_unordered(get_all_base_probabilities_matrices_one_fragment, arguments_list):
        out_dict[dict_element['name']] = dict_element
        i += 1
        if i % how_often_print == 0:
            toc = time.time()
            print("Processed %d fragments. Time spent: %d" % (i, toc-tic))
    return out_dict


def caluclate_base_pairing_differences(matrices_dict_element):
    out_dict = {}

    out_dict['sum_unconstr'] = np.sum(matrices_dict_element['unconstr'])

    out_dict['dms_1_unconstr'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS1_unconstr'], matrices_dict_element['unconstr'])
    out_dict['dms_2_unconstr'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS2_unconstr'], matrices_dict_element['unconstr'])

    out_dict['major_loop_unconstr'] = difference_normalized_by_second_matrix(
        matrices_dict_element['major_loop'], matrices_dict_element['unconstr'])
    out_dict['second_loop_unconstr'] = difference_normalized_by_second_matrix(
        matrices_dict_element['second_loop'], matrices_dict_element['unconstr'])

    out_dict['dms_1_major_loop'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS1_major_loop'], matrices_dict_element['major_loop'])
    out_dict['dms_1_second_loop'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS1_second_loop'], matrices_dict_element['second_loop'])
    out_dict['dms_2_major_loop'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS2_major_loop'], matrices_dict_element['major_loop'])
    out_dict['dms_2_second_loop'] = difference_normalized_by_second_matrix(
        matrices_dict_element['DMS2_second_loop'], matrices_dict_element['second_loop'])
    return out_dict


def discretize_count_dms_profile(dms_profile,
                               max_value,
                               n_bins):
    bins_borders = np.zeros(n_bins + 1, dtype = np.float)
    bins_borders[:-1] = np.linspace(0, max_value, n_bins)
    bins_borders[-1] = 1000
    dms_profile_non_neg = dms_profile.copy()
    dms_profile_non_neg[dms_profile_non_neg < 0] = 0
    histogram = np.histogram(dms_profile_non_neg, bins = bins_borders)
    return histogram


# copied from pyteiser MI
def entropy_empirical(counts, total_number, base=None):
    probs = np.divide(counts, total_number)
    ent = 0.
    base = math.e if base is None else base

    for i in probs:
        if i == 0: # np.isclose is not supported by numba
            continue
        ent -= i * math.log(i) / math.log(base)

    return ent


def profile_entropy(dms_profile,
                    max_value = 1.5,
                    n_bins = 16
                    ):
    histogram, bins = discretize_count_dms_profile(
                           dms_profile,
                           max_value = max_value,
                           n_bins = n_bins)
    entropy = entropy_empirical(histogram, histogram.sum())
    return entropy

def write_common_unwound_structure(fragment_name, temp_files_folder, ct_filenames):
    new_ct_file = os.path.join(temp_files_folder, fragment_name + '_basepairs_common_between_two_major_consensuses.ct')
    assert len(ct_filenames) == 2, "Too many .ct files! Error\n\n\n"
    with open(new_ct_file, 'w') as wf:
        with open(ct_filenames[0], 'r') as rf1:
            with open(ct_filenames[1], 'r') as rf2:
                for line1, line2 in zip(rf1, rf2):
                    if 'ENERGY' in line1 and 'ENERGY' in line2:
                        wf.write(line1)
                        continue
                    if line1 == line2:
                        wf.write(line1)
                    else:
                        line_single_spaces = re.sub(' +', ' ', line1)
                        line_splitted = line_single_spaces.split(' ')
                        line_splitted[4] = '0'
                        str_to_write = ''
                        str_to_write += ' ' * (5 - len(line_splitted[1]))
                        str_to_write += line_splitted[1]
                        str_to_write += ' '+ line_splitted[2]
                        str_to_write += ' ' * (8 - len(line_splitted[3]))
                        str_to_write += line_splitted[3]
                        str_to_write += ' ' * (5 - len(line_splitted[4]))
                        str_to_write += line_splitted[4]
                        str_to_write += '    0'
                        str_to_write += ' ' * (6 - len(line_splitted[6]))
                        str_to_write += line_splitted[6]
                        wf.write(str_to_write)
    return new_ct_file


def ctfile_to_dotbracket(inp_filename, RNAstructure_path):
    out_filename = inp_filename.replace('.ct','.txt')
    ct2dot_path = os.path.join(RNAstructure_path, "ct2dot")
    conversion_arguments = [ct2dot_path, inp_filename, '1', out_filename]
    ct_result = run(args=conversion_arguments, stdout=PIPE, stderr=PIPE)
    # print("ct_result", ct_result.stderr)
    with open(out_filename, 'r') as rf:
        splitfile = rf.read().split('\n')
        bracket_notation = splitfile[2]
    return bracket_notation, out_filename


def calculate_free_energies(infiles_names, RNAstructure_path):
    efn2_path = os.path.join(RNAstructure_path, "efn2")
    outfiles_names = [x.replace('.ct', '_free_energy.ct') for x in infiles_names]
    energies_calculated = []

    for fin, fout in zip(infiles_names, outfiles_names):
        calculation_arguments = [efn2_path, fin, fout]
        calc_result = run(args=calculation_arguments, stdout=PIPE, stderr=PIPE)
        # print("calc", calc_result.stderr)
        with open(fout, 'r') as rf:
            current_energy = float(rf.read().split(' ')[-1])
            energies_calculated.append(current_energy)

    free_energies_dict = {'major_loop': energies_calculated[0],
                          'second_loop': energies_calculated[1],
                          'saddle': energies_calculated[2]}

    return free_energies_dict, outfiles_names


def fold_perturbation_shape_MEA(fasta_file,
                                shape_file,
                                temp_folder,
                                RNAstructure_path):
    RNAstructure_parameters_data_path = os.path.join(RNAstructure_path, 'data_tables')
    RNAstructure_path = os.path.join(RNAstructure_path, 'exe')
    # RNAstructure won't work without DATAPATH environmental variable set to the right folder
    # I can't set it with subprocess.Popen for reasons described here: https://stackoverflow.com/questions/3929319/subprocess-module-errors-with-export-in-python-on-linux
    os.putenv("DATAPATH", RNAstructure_parameters_data_path)

    partition_command = os.path.join(RNAstructure_path, "partition")
    MEA_command = os.path.join(RNAstructure_path, "MaxExpect")

    fasta_sample_name = os.path.basename(fasta_file).replace('.fa','')
    shape_sample_name = os.path.basename(shape_file).replace('_shape.dat', '')
    full_sample_name = "%s_%s" % (fasta_sample_name, shape_sample_name)

    pfs_filename = os.path.join(temp_folder, "%s_partition.pfs" % full_sample_name)
    ct_filename = os.path.join(temp_folder, "%s_MEA.ct" % full_sample_name)
    filenames_to_remove = [pfs_filename, ct_filename]

    partition_arguments = [partition_command, fasta_file, pfs_filename]
    MEA_arguments = [MEA_command, pfs_filename, ct_filename]

    if not shape_file is None:
        partition_arguments.append('--SHAPE')
        partition_arguments.append(shape_file)

    partition_result = run(args=partition_arguments, stdout=PIPE, stderr=PIPE, cwd=temp_folder)
    MEA_result = run(args=MEA_arguments, stdout=PIPE, stderr=PIPE, cwd=temp_folder)
    MEA_bracket_notation, dot_filename = ctfile_to_dotbracket(ct_filename, RNAstructure_path)
    filenames_to_remove.append(dot_filename)
    string_to_write = ">%s\n%s\n" % (full_sample_name, MEA_bracket_notation)
    for fn in sorted(list(set(filenames_to_remove))):
        os.remove(fn)
    return string_to_write


def fold_perturbation_shape_MFE(fasta_file,
                                shape_file,
                                temp_folder,
                                RNAstructure_path):
    RNAstructure_parameters_data_path = os.path.join(RNAstructure_path, 'data_tables')
    RNAstructure_path = os.path.join(RNAstructure_path, 'exe')
    # RNAstructure won't work without DATAPATH environmental variable set to the right folder
    # I can't set it with subprocess.Popen for reasons described here: https://stackoverflow.com/questions/3929319/subprocess-module-errors-with-export-in-python-on-linux
    os.putenv("DATAPATH", RNAstructure_parameters_data_path)

    fold_command = os.path.join(RNAstructure_path, "Fold")

    fasta_sample_name = os.path.basename(fasta_file).replace('.fa','')
    shape_sample_name = os.path.basename(shape_file).replace('_shape.dat', '')
    full_sample_name = "%s_%s" % (fasta_sample_name, shape_sample_name)

    ct_filename = os.path.join(temp_folder, "%s_MEA.ct" % full_sample_name)
    filenames_to_remove = [ct_filename]

    fold_arguments = [fold_command, fasta_file, ct_filename]

    if not shape_file is None:
        fold_arguments.append('--SHAPE')
        fold_arguments.append(shape_file)

    fold_result = run(args=fold_arguments, stdout=PIPE, stderr=PIPE, cwd=temp_folder)
    MEA_bracket_notation, dot_filename = ctfile_to_dotbracket(ct_filename, RNAstructure_path)
    filenames_to_remove.append(dot_filename)
    string_to_write = ">%s\n%s\n" % (full_sample_name, MEA_bracket_notation)
    for fn in sorted(list(set(filenames_to_remove))):
        os.remove(fn)
    return string_to_write


def conflicting_loops_scores_output_parser(inp_scores_file_loc):
    init_pars_dict_loc = {}
    with open(inp_scores_file_loc, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n$$$\n')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            # print(entry)
            splitted_entry = entry.split('\n')
            if splitted_entry[1] == 'ERROR':
                continue
            if len(splitted_entry) <= 2:
                continue
            fragment_name = splitted_entry[0]
            basepairs_number = re.findall("\d+", splitted_entry[1])[0]
            basepairs_number = float(basepairs_number)
            init_pars_dict_loc[fragment_name] = {}
            init_pars_dict_loc[fragment_name]['basepairs_number'] = basepairs_number
            init_pars_dict_loc[fragment_name]['pairs_loops'] = {}
            format_counter = -1
            for line in splitted_entry[2:]:
                if line.startswith('Overlapping part of loops'):
                    pair_number = re.findall("\d+", line)[0]
                    pair_number = int(pair_number)
                    init_pars_dict_loc[fragment_name]['pairs_loops'][pair_number] = {}
                    format_counter = 0
                    continue

                current_loop = 0
                if format_counter == 0:
                    current_loop = 0
                    format_counter = 1
                elif format_counter == 1:
                    current_loop = 1
                    format_counter = 2
                else:
                    print("Error!")
                    break

                numbers_in_string = re.findall(r'\d+\.*\d*', line)
                numbers_in_string = [float(x) for x in numbers_in_string]
                init_pars_dict_loc[fragment_name]['pairs_loops'][pair_number][current_loop] = \
                    {"forward": (numbers_in_string[0], numbers_in_string[1]), \
                     "backward": (numbers_in_string[2], numbers_in_string[3]),
                     "average_MI": numbers_in_string[4]}
    return init_pars_dict_loc


def get_all_the_fragments_from_fasta_file(input_fasta_file_loc):
    fragments_seq_dict_loc = {}
    with open(input_fasta_file_loc, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n>')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            index = entry.find('\n')
            fragment_name_loc = entry[0:index]
            fragment_name_loc = fragment_name_loc.replace('>', '')
            sequence = entry[index:].replace('\n', '')
            fragments_seq_dict_loc[fragment_name_loc] = sequence
    return fragments_seq_dict_loc


def get_two_major_conflicting_loops(dict_element):
    major_loop, major_loop_MI = get_major_loop(dict_element)
    second_loop = get_max_loop_conflicting_with_the_major_one(dict_element, major_loop, major_loop_MI)
    loops_dict_fr = {'major':major_loop, 'second':second_loop}
    return loops_dict_fr


def get_major_loop(init_pars_dict_entry):
    major_loop = {}
    major_loop_MI = 0
    for num1 in init_pars_dict_entry['pairs_loops']:
        for num2 in init_pars_dict_entry['pairs_loops'][num1]:
            current_loop = init_pars_dict_entry['pairs_loops'][num1][num2]
            current_loop_MI = current_loop['average_MI']
            if current_loop_MI > major_loop_MI:
                major_loop_MI = current_loop_MI
                major_loop = current_loop
    return major_loop, major_loop_MI


def get_max_loop_conflicting_with_the_major_one(init_pars_dict_entry, major_loop, major_loop_MI):
    second_loop = {}
    second_loop_MI = 0
    major_loop_characteristic_string = get_characteristic_string(major_loop)

    for num1 in init_pars_dict_entry['pairs_loops']:
        characteristic_strings_dict = {}
        for num2 in init_pars_dict_entry['pairs_loops'][num1]:
            current_loop = init_pars_dict_entry['pairs_loops'][num1][num2]
            characteristic_string = get_characteristic_string(current_loop)
            characteristic_strings_dict[characteristic_string] = current_loop
        if major_loop_characteristic_string in characteristic_strings_dict:
            del characteristic_strings_dict[major_loop_characteristic_string]
            keys_list = sorted(list(characteristic_strings_dict.keys()))
            current_loop = characteristic_strings_dict[keys_list[0]]
            current_loop_MI = current_loop['average_MI']
            if len(keys_list) != 1:
                print("Error!")
            if current_loop_MI > second_loop_MI:
                second_loop = current_loop
                second_loop_MI = current_loop_MI

    return second_loop


def make_sting_to_write(fragment_name, fragment_sequence, free_energies_dict, consensus_structures_dict):
    string_to_write = '%s\n%s\n' % (fragment_name, fragment_sequence)
    string_to_write += 'major_loop\n%s\nsecond_loop\n%s\ncommon\n%s\n' % (consensus_structures_dict['major'], consensus_structures_dict['second'], consensus_structures_dict['common'])
    string_to_write += 'major_loop: %f, second_loop: %f, saddle: %f\n' % (free_energies_dict['major_loop'], free_energies_dict['second_loop'], free_energies_dict['saddle'])
    string_to_write += '$$$\n'
    return string_to_write
