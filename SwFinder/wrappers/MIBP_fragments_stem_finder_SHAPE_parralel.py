import argparse
import os
import sys
import subprocess
import numpy as np
import re
from datetime import datetime
import multiprocessing
import pickle
import random
import string


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="the fasta file with the target sequences", type=str)
    parser.add_argument("--shape_dict", help="dictionary with SHAPE/DMS values", type=str)
    parser.add_argument("--temp_files_folder", help="for files like RNA_name*MI.txt", type=str)
    parser.add_argument("-o", help="output folder", type=str)
    parser.add_argument("--eta", help="throw away pairs with prob less than eta or greater than 1-eta when calculating MIBP (see luan lin's thesis for why)", type=float)
    parser.add_argument("--delta", help="minimum bp prob", type=float)
    parser.add_argument("--MI_max_ratio_filtering_constant", help="what portion of maximal MI should a pair carry to not be thrown away", type=float)
    parser.add_argument("--min_stretch_length",help="", type=int)
    parser.add_argument("--min_stretch_overlap", help="", type=int)
    parser.add_argument("--RNAstructure_path", help="path to RNAstructure program", type=str)
    parser.add_argument("--num_samples", help="NUMBER OF SAMPLES IN STOCHASTIC SAMPLING", type=int)
    parser.add_argument("--num_processes", help="number of processes to run in parralel", type=int)
    parser.add_argument("--how_often_print", help="print when each N-th fragment is done", type=int)


    parser.set_defaults(
                        shape_dict = None,
                        temp_files_folder='/avicenna/khorms/projects/SNIP_switchers/MIBP_python/temp_files',
                        RNAstructure_path='/avicenna/khorms/programs/RNAstructure',
                        eta = 0.01,
                        delta = 0.001,
                        num_samples = 1000,
                        num_processes = 20,
                        min_stretch_length=3,
                        min_stretch_overlap=3,
                        MI_max_ratio_filtering_constant=0.5,
                        how_often_print=10
                        )
    args = parser.parse_args()
    return args


def define_variables():
    args = handler()

    global RT
    global NUM_SAMPLES
    global NUM_PROCESSES
    global ETA
    global DELTA
    global MIN_STRETCH_LENGTH
    global MIN_STRETCH_OVERLAP
    global MI_MAX_RATIO_FILTERING_CONSTANT
    global HOW_OFTEN_TO_PRINT

    RT = 0.61633
    ETA = args.eta
    DELTA = args.delta
    NUM_SAMPLES = args.num_samples
    NUM_PROCESSES = args.num_processes
    MIN_STRETCH_LENGTH = args.min_stretch_length
    MIN_STRETCH_OVERLAP = args.min_stretch_overlap
    MI_MAX_RATIO_FILTERING_CONSTANT = args.MI_max_ratio_filtering_constant
    HOW_OFTEN_TO_PRINT = args.how_often_print

    input_fasta_file = args.f
    short_filename = input_fasta_file.split('/')[-1].split('.')[0]
    input_shape_dict = pickle.load(open(args.shape_dict, 'rb'))
    output_folder = args.o
    output_filename = os.path.join(output_folder, short_filename + '_parralel_MIBP_output.txt')
    output_timing_filename = os.path.join(output_folder, short_filename + '_parralel_MIBP_timing.txt')

    global temp_files_folder
    global RNAstructure_parameters_data_path
    global partition_program_path
    global stochastic_program_path

    temp_files_folder = args.temp_files_folder
    RNAstructure_path = args.RNAstructure_path

    RNAstructure_parameters_data_path = os.path.join(RNAstructure_path, 'data_tables')
    partition_program_path = os.path.join(RNAstructure_path, 'exe/partition')
    stochastic_program_path = os.path.join(RNAstructure_path, 'exe/stochastic')


    return input_fasta_file, input_shape_dict, output_filename, output_timing_filename


def run_RNAstructure_precalculations(fragment_sequence_file, filenames_dict):
    confile = filenames_dict['confile']
    create_constraint_file([], [], [], [], confile)
    pfsfile = filenames_dict['pfsfile']
    ctfile = filenames_dict['ctfile']
    shapefile = filenames_dict['shape']

    ## Calculate partition function (Q) of the entire ensemble:
    # (Q will be used to calculate the probability mass of the "sub-spaces")

    datapath_command = "export DATAPATH=%s ;" % (RNAstructure_parameters_data_path)

    partition_command = "%s -c %s %s %s" % (partition_program_path, confile, fragment_sequence_file, pfsfile)
    if not shapefile is None:
        partition_command += ' --SHAPE %s' % (shapefile)
    partition_command = datapath_command + partition_command
    subprocess.call(partition_command, shell=True)

    stochastic_command = "%s %s %s" % (stochastic_program_path, pfsfile, ctfile)
    stochastic_command = datapath_command + stochastic_command
    subprocess.call(stochastic_command, shell=True)


def energy2partition(inp_energy):
    # Converts from "ensemble energy" (see
    # http://rna.urmc.rochester.edu/Text/EnsembleEnergy.html) to partition
    # function, Q
    return np.exp(inp_energy/(-1*RT))


def make_MI_table(fragment_name_loc, filenames_dict, sequence):
    ctfile = filenames_dict['ctfile']
    structures = readct(ctfile, len(sequence))

    # get mutual information
    max_mi, basepairs_sum_MIs_raw, base_pairs = rna_mi3(structures)

    the_pairs_we_consider = basepairs_sum_MIs_raw[:,2] > MI_MAX_RATIO_FILTERING_CONSTANT * max_mi
    basepairs_sum_MIs_filtered = np.copy(basepairs_sum_MIs_raw)
    basepairs_sum_MIs_filtered = basepairs_sum_MIs_filtered[the_pairs_we_consider, :]

    stems_dict = look_for_stretches(basepairs_sum_MIs_filtered)
    average_stretch_mi_dict = average_stems_MI(stems_dict)
    number_of_prob_basepairs = base_pairs.shape[0]
    stems_boundaries = find_stem_boundaries(stems_dict)
    overlapping_stems = find_overlapping_stems(stems_dict, stems_boundaries)
    string_to_write = write_overlapping_stems_to_file(stems_dict, overlapping_stems,
                                                      stems_boundaries, average_stretch_mi_dict,
                                                    number_of_prob_basepairs, fragment_name_loc)
    return string_to_write


def look_for_stretches(inp_basepairs_table):
    stems_dict = {}
    max_index_value = inp_basepairs_table.shape[0]
    marked_ones = np.zeros(inp_basepairs_table.shape[0], dtype=bool)
    for index in range(max_index_value):
        if marked_ones[index]:
            continue
        stems_dict[index] = 0
        base = inp_basepairs_table[index, 0]
        pair = inp_basepairs_table[index, 1]
        current_array = inp_basepairs_table[index,:]
        sub_bool_final = np.copy(marked_ones)
        for sub_index in range(index+1, max_index_value):
            sub_table = inp_basepairs_table[inp_basepairs_table[:,0] == base + sub_index - index, :]
            sub_bool = inp_basepairs_table[:,0] == base + sub_index - index
            if sub_table.shape[0] == 0:
                break
            sub_table_2 = sub_table[sub_table[:,1] == pair - sub_index + index, :]
            sub_bool_2 = inp_basepairs_table[:,1] == pair - sub_index + index
            sub_bool_combined = sub_bool & sub_bool_2
            sub_bool_final = sub_bool_combined | sub_bool_final
            if sub_table_2.shape[0] == 0:
                break
            current_array = np.vstack((current_array, sub_table_2))
            stems_dict[index] += 1
        if stems_dict[index] < MIN_STRETCH_LENGTH:
            del stems_dict[index]
        else:
            stems_dict[index] = current_array
            marked_ones = marked_ones | sub_bool_final
            #print(current_array)

    return stems_dict


def average_stems_MI(stems_dict):
    average_mi_dict = {}
    for stretch in stems_dict:
        mean_MI = np.mean(stems_dict[stretch][:,2])
        # number_of_prob_basepairs = base_pairs.shape[0]
        # MI_to_bp_number_value = round(mean_MI / np.log2(number_of_prob_basepairs), 2)
        average_mi_dict[stretch] = mean_MI
    return average_mi_dict


def find_stem_boundaries(stems_dict):
    stems_boundaries = {}
    for stretch in stems_dict:
        stretch_forward = (int(stems_dict[stretch][0,0]), int(stems_dict[stretch][-1,0]))
        stretch_backward = (int(stems_dict[stretch][-1, 1]), int(stems_dict[stretch][0, 1]))
        stems_boundaries[stretch] = {"forward":stretch_forward, "backward":stretch_backward}
    return stems_boundaries


def find_overlapping_stems(stems_dict, stems_boundaries):
    stretch_indices = sorted(list(stems_dict.keys()))
    overlapping_stems = np.zeros((len(stretch_indices), len(stretch_indices)))

    for p in range(len(stretch_indices)):
        for q in range(len(stretch_indices)):
            i = stretch_indices[p]
            k = stretch_indices[q]
            if i == k:
                continue
            forward_forward_overlap = get_overlap(stems_boundaries[i]["forward"], stems_boundaries[k]["forward"])
            backward_backward_overlap = get_overlap(stems_boundaries[i]["backward"], stems_boundaries[k]["backward"])
            forward_backward_overlap = get_overlap(stems_boundaries[i]["forward"], stems_boundaries[k]["backward"])
            backward_forward_overlap = get_overlap(stems_boundaries[i]["backward"], stems_boundaries[k]["forward"])

            if (forward_forward_overlap >= MIN_STRETCH_OVERLAP) or \
                    (backward_backward_overlap >= MIN_STRETCH_OVERLAP) or \
                    (forward_backward_overlap >= MIN_STRETCH_OVERLAP) or \
                    (backward_forward_overlap >= MIN_STRETCH_OVERLAP):
                overlapping_stems[p, q] = 1
            else:
                overlapping_stems[p, q] = 0

    overlapping_stems = np.triu(overlapping_stems)
    return overlapping_stems


def write_overlapping_stems_to_file(stems_dict, overlapping_stems, stems_boundaries,
                                    average_stretch_mi_dict, number_of_prob_basepairs, fragment_name_loc):

    stretch_indices = sorted(list(stems_dict.keys()))
    string_to_write = ''
    string_to_write += '%s\n' % (fragment_name_loc)
    string_to_write += "Number of possible basepairs: %d\n" % (number_of_prob_basepairs)
    counter = 0
    for p in range(len(stretch_indices)):
        for q in range(len(stretch_indices)):
            i = stretch_indices[p]
            k = stretch_indices[q]
            if overlapping_stems[p,q] == 1:
                counter += 1
                string_to_write += "Overlapping part of loops %d:\n" % (counter)
                first_stem_string = str(stems_boundaries[i]["forward"]) + ' | ' + str(stems_boundaries[i]["backward"])
                first_stem_string += ' | Average MI: %f\n' % (average_stretch_mi_dict[i])
                string_to_write += first_stem_string
                second_stem_string = str(stems_boundaries[k]["forward"]) + ' | ' + str(stems_boundaries[k]["backward"])
                second_stem_string += ' | Average MI: %f\n' % (average_stretch_mi_dict[k])
                string_to_write += second_stem_string

    string_to_write += '$$$\n'
    return string_to_write


def rna_mi3(samples):
    # This function finds the mutual information (MI) of all pairs of nucleotides
    # from a sampled collection of secondary structures.
    # INPUT: structures, a matrix of base pairs (see readct.m)
    # OUTPUTS: max_mi, a num, the maximum average MI of any pair of nucleotides
    #            with all other pairs of nucleotides.
    #          max_pair, a length-2 array of ints, the pair of nucleotides
    #            achieving the max average MI
    #          avg_mi, an m x 3 matrix, where the first two entries of a row
    #            are a pair of nucleotides, and the third entry is the average
    #              MI this pair has with all other pairs.

    counts = count_from_struct(samples)  # make a table of base pairs for each structure

    base_pairs = pair_probs2(counts)  # get the base pair probs
    base_pairs = base_pairs[base_pairs[:, 2] > DELTA,:]  # filter out pairs with low prob
    # (see Luan Lin's thesis for theoretical justification of filtering out low
    # probs)

    # do the real work here
    # loop through all pairs of base pairs
    # this can probably be made faster
    mi = np.zeros((len(base_pairs), len(base_pairs)))
    for i in range(base_pairs.shape[0] - 1):
        if (base_pairs[i, 2] > ETA) and (base_pairs[i, 2] < 1 - ETA):  # get rid of low entropy pairs
            p1_count = (samples[:, [int(x) for x in base_pairs[i, 0].flat]] == base_pairs[i, 1]) # get a vector of all structures that have this bp
            for j in range((i+1),base_pairs.shape[0]):
                if (base_pairs[j, 2] > ETA) and (base_pairs[j, 2] < 1 - ETA):
                    p2_count = (samples[:, [int(x) for x in base_pairs[j, 0].flat]] == base_pairs[j, 1])
                    mi[i, j] = mutual_info(p1_count, p2_count)  # get the MI between these two bps using the structures
                    mi[j, i] = mi[i, j]  # MI is symmetric

    # get the max average MI, throwing away pairs with low entropy
    # for why we throw away pairs with low entropy, see Luan Lin's thesis
    c = 1
    sum_mi = np.zeros((base_pairs.shape[0], 3))
    for i in range(len(base_pairs)):
        if (base_pairs[i, 2] < 1 - ETA) and (base_pairs[i, 2] > ETA):
            sum_mi[c, 0:2] = base_pairs[i, 0:2]
            sum_mi[c, 2] = sum(mi[i,:])
            c += 1

    boolean_vector_to_filter = sum_mi[:, 2] != 0
    sum_mi = sum_mi[boolean_vector_to_filter, :] # filter out the ones that have zero in the second column? khorms sum_mi(sum_mi(:,3) == 0, :) = [];
    # mi = mi[boolean_vector_to_filter, :] # filter them from this table too
    # mi = np.transpose(mi)[boolean_vector_to_filter, :]
    # mi = np.transpose(mi)

    max_mi = max(sum_mi[:, 2])

    return max_mi, sum_mi, base_pairs # return maximal MI value, table of all sums of MI and a table with all probable basepairs


def mutual_info(p1, p2):
    # calculate MI of two binary random variables from samples
    # INPUTS: An array that gives the value of the first random variable in
    # each sample, an array that gives the value of the second random variable in
    # each sample.
    # Count how many time each r.v. is 1
    p1_count = np.sum(p1)
    p2_count = np.sum(p2)

    # Marginal for first r.v.
    prob1 = np.array([p1_count/float(NUM_SAMPLES), (NUM_SAMPLES - p1_count)/float(NUM_SAMPLES)])

    # Marginal for second r.v.
    prob2 = np.array([p2_count/float(NUM_SAMPLES), (NUM_SAMPLES - p2_count)/float(NUM_SAMPLES)])

    # Joint distribution
    #SUM(condition) doesn't work in python!
    prob3 = np.array([np.sum(np.logical_and(p1 == 0, p2 == 0))/float(NUM_SAMPLES), np.sum(np.logical_and(p1 == 0, p2 == 1))/float(NUM_SAMPLES),
                      np.sum(np.logical_and(p1 == 1, p2 == 0))/float(NUM_SAMPLES), np.sum(np.logical_and(p1 == 1, p2 == 1))/float(NUM_SAMPLES)])

    # I(X;Y) = H(X) + H(Y) - H(X,Y)

    mi = shannon_entropy(prob1) + shannon_entropy(prob2) - shannon_entropy(prob3)
    return mi


def shannon_entropy(x):
    # This function computes the entropy of one or more discrete random variables. It is
    # used to find the entropy of all base pairs.
    # INPUT: x, a k x m matrix, where k is the number of values the variables
    #          can on (so k = 2 in the case of base pairs), and m is the number
    #          of variables. Each entry in a column is the probability that the
    #          variable takes on the corresponding value.
    # OUTPUT: entropy, an m-length array, where the i'th entry is the entropy
    #           of the i'th random variable.

    # mask = np.ones(x.shape, dtype=bool) # a hack to get around NaN
    # mask[np.flatnonzero(x)] = False
    # x[mask] = 1
    # entropy = -sum(np.multiply(x, np.log2(x)))
    # return entropy

    log_x = np.log2(x)
    log_x[np.isneginf(log_x)] = 0
    entropy = -np.sum(np.multiply(x, log_x))
    return entropy


def pair_probs2(counts):
    # Converts the counts returned by count_pairs.m into base pair
    # probabilities.
    # INPUT: structure, an L x L array (see readct.m) recording base
    #          pairs
    # OUTPUT: base_pairs, an m x 3 matrix, where m is the number of distinct base pairs
    #           that were observed in the sampled structures. base_pairs(j, 1)
    #           and base_pairs(j, 2) are the nucleotides of the j'th pair of
    #           bases. nucleotides(j, 3) is the sample probability of this base
    #           pair forming

    r = np.where(counts > 0)[0]
    c = np.where(counts > 0)[1]
    base_pairs = np.zeros((len(r),3))
    base_pairs[:,0] = r
    base_pairs[:,1] = c
    base_pairs[:,2] = counts[counts>0]/float(NUM_SAMPLES)
    if base_pairs[0,0] > base_pairs[1,0]:
        store_vect = base_pairs[:,0]
        base_pairs[:, 0] = base_pairs[:,1]
        base_pairs[:, 1] = store_vect

    return base_pairs


def count_from_struct(samples):
    # Count how many times each base pair appears in a sample.
    # INPUT: A 1000 x RNA_LENGTH array describing which pairs apear in which
    # samples
    # OUTPUT: An RNA_LENGTH x RNA_LENGTH array with base pair counts.
    counts = np.zeros((samples.shape[1], samples.shape[1]))
    for i in range(samples.shape[0]): # Loop through samples
        for j in range(samples.shape[1]): # Loop through sequence
            pairedto = samples[i,j]
            if j < pairedto:
                # Add to the array if there is a base pair.
                counts[j,int(pairedto)] += 1
    return counts


def create_constraint_file(forced_paired, forced_unpaired, forced_pairs, forced_nonpairs, file_name_inp):
    with open(file_name_inp, 'w') as wf:
        wf.write('DS:\n')
        for i in range(0, len(forced_paired)):
            wf.write(str(int(forced_paired[i] + 1)) + '\n')
        wf.write('-1\nSS:\n')
        for i in range(0, len(forced_unpaired)):
            wf.write(str(int(forced_unpaired[i] + 1)) + '\n')
        wf.write('-1\nMod:\n-1\nPairs:\n')
        for i in range(0, len(forced_pairs)):
            wf.write(str(int(forced_pairs[i][0] + 1)) + ' ' + str(int(forced_pairs[i][1] + 1)) + '\n')
        wf.write('-1 -1\nFMN:\n-1\nForbids:\n')
        for i in range(0, len(forced_nonpairs)):
            wf.write(str(int(forced_nonpairs[i][0] + 1)) + ' ' + str(int(forced_nonpairs[i][1] + 1)) + '\n')
        wf.write('-1 -1')


def readct(ctfile, seq_length):
    # Get column 5 from the ct file. Output is a NUM_SAMPLES x seq_length
    # array samples(i,j) is pair of base j in sample i or 0 if it is
    # unpaired.

    samples = np.zeros((NUM_SAMPLES,seq_length))

    with open(ctfile, 'r') as rf:
        counter = 0
        num_samples_counter = 0
        for line in rf:
            counter += 1
            if counter % (seq_length + 1) == 1:
                continue
            if counter > (NUM_SAMPLES * (seq_length+1)):
                break
            if counter % (seq_length + 1) == 0:
                num_samples_counter += 1
            line_single_spaces = re.sub(' +',' ',line)
            line_splitted = line_single_spaces.split(' ')
            pos1 = int(line_splitted[1]) - 1
            pos2 = int(line_splitted[5])
            # if pos2 > 0:
            #     pos2 -= 1

            pos2 -= 1

            samples[num_samples_counter - 1,pos1] = pos2
    return samples


def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def write_fragment_sequence_file(fragment_name_loc, seq):
    fragment_sequence_file = os.path.join(temp_files_folder, fragment_name_loc + '_sequence.fa')
    with open(fragment_sequence_file, 'w') as wf:
        wf.write('>%s\n%s\n' % (fragment_name_loc, seq))
    return fragment_sequence_file


def write_shape_constraints(filename, shape_array):
    strings_to_write_list = []
    for i in range(shape_array.shape[0]):
        k = i + 1 # move to 1-based coordinate system
        strings_to_write_list.append("%d\t%.5f" % (k, shape_array[i]))

    string_to_write = "\n".join(strings_to_write_list)
    with open(filename, 'w') as wf:
        wf.write(string_to_write)


def create_files_required(fragment_name_loc, shape_profile):
    rand_string = ''.join(random.choices(string.ascii_uppercase + string.digits, k=7))
    filenames_dict = {}
    filenames_dict['confile'] = os.path.join(temp_files_folder, '%s_%s_.CON' % (fragment_name_loc, rand_string))
    filenames_dict['pfsfile'] = os.path.join(temp_files_folder, '%s_%s_.pfs' % (fragment_name_loc, rand_string))
    filenames_dict['ctfile'] = os.path.join(temp_files_folder, '%s_%s_.ct' % (fragment_name_loc, rand_string))

    for fn_var in filenames_dict:
        with open(filenames_dict[fn_var], 'w') as wf:
            pass

    if not shape_profile is None:
        filenames_dict['shape'] = os.path.join(temp_files_folder, '%s_%s_shape.dot' % (fragment_name_loc, rand_string))
        write_shape_constraints(filenames_dict['shape'], shape_profile)
    else:
        filenames_dict['shape'] = None

    return filenames_dict


def delete_files_required(filenames_dict, fragment_name_loc):
    for fn_var in filenames_dict:
        os.remove(filenames_dict[fn_var])
    fragment_sequence_file = os.path.join(temp_files_folder, fragment_name_loc + '_sequence.fa')
    os.remove(fragment_sequence_file)


def worker_for_parallel_implementation(curr_tuple):
    fragment_name_loc = curr_tuple[0]
    sequence = curr_tuple[1]
    if len(curr_tuple) > 2:
        shape_profile = curr_tuple[2]
    else:
        shape_profile = None
    fragment_sequence_file = write_fragment_sequence_file(fragment_name_loc, sequence)
    filenames_dict = create_files_required(fragment_name_loc, shape_profile)
    try:
        run_RNAstructure_precalculations(fragment_sequence_file, filenames_dict)
        string_to_write = make_MI_table(fragment_name_loc, filenames_dict, sequence)
    except:
        return '%s\nERROR\n$$$\n' % (fragment_name_loc)
    delete_files_required(filenames_dict, fragment_name_loc)
    return string_to_write


def make_an_iterable_for_multiprocessing(input_fasta_file, input_shape_dict):
    iterator_list_of_inputs = []
    with open(input_fasta_file, 'r') as rf:
        bigline = rf.read()
        bigarray = bigline.split('\n>')
        bigarray = [x for x in bigarray if x != '']
        for entry in bigarray:
            index = entry.find('\n')
            fragment_name_loc = entry[0:index]
            fragment_name_loc = fragment_name_loc.replace('>', '')
            sequence = entry[index:].replace('\n', '')
            curr_annotated_list = [fragment_name_loc, sequence]

            if not input_shape_dict is None:
                if fragment_name_loc not in input_shape_dict:
                    print("Fragment %s not in SHAPE profiles dict, even though SHAPE values have been provided" %
                          fragment_name_loc)
                    sys.exit(1)
                curr_shape_profile = input_shape_dict[fragment_name_loc]
                assert len(sequence) == curr_shape_profile.shape[0]
                curr_annotated_list.append(curr_shape_profile)

            curr_tuple = tuple(curr_annotated_list)
            iterator_list_of_inputs.append(curr_tuple)

    return iterator_list_of_inputs


def mp_handler(input_fasta_file, input_shape_dict, output_filename, output_timing_filename):
    pool = multiprocessing.Pool(NUM_PROCESSES)
    iterator_list_of_inputs = make_an_iterable_for_multiprocessing(input_fasta_file, input_shape_dict)
    counter = 0
    with open(output_filename, 'w') as wf:
        with open(output_timing_filename, 'w') as write_timing:
            print("Started")
            write_timing.write("Beginning the calculation:\t" +
                               datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
            for string_to_write in pool.imap_unordered(worker_for_parallel_implementation, iterator_list_of_inputs):
                print('Done: fragment number ', counter)
                counter += 1
                wf.write(string_to_write)
                wf.flush()
                if counter % HOW_OFTEN_TO_PRINT == 0:
                    write_timing.write(("%d fragments have been processed:\t" % (counter)) +
                                       datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\n')
                    write_timing.flush()


def main():
    input_fasta_file, input_shape_dict, output_filename, output_timing_filename = define_variables()
    mp_handler(input_fasta_file, input_shape_dict, output_filename, output_timing_filename)


if __name__ == "__main__":
    main()