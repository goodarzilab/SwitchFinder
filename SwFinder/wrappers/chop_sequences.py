import argparse
import IO as IO
import numpy as np

def handler(raw_args = None):
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", help="the fasta file with the target sequences", type=str)
    parser.add_argument("-o", help="output filename", type=str)
    parser.add_argument("--length", help="fragment length", type=int)


    parser.set_defaults(
                        f = "",
                        o = "",
                        length = 186
                        )
    args = parser.parse_args(raw_args)
    return args


def make_intervals_for_rw_separation_v2(utr_length, length_limit):
    if utr_length > length_limit:
        half_length_limit = length_limit // 2
        number_of_fragments = utr_length // half_length_limit
        leftover_portion = utr_length % length_limit

        if leftover_portion == 0:
            cut_coordinates = np.zeros((number_of_fragments - 1, 3), dtype=int)
        else:
            cut_coordinates = np.zeros((number_of_fragments, 3), dtype=int)

        for k in range(0, number_of_fragments - 1):
            start_coordinate = k * length_limit // 2
            end_coordinate = (k + 2) * length_limit // 2
            cut_coordinates[k, :] = [k + 1, start_coordinate, end_coordinate]

        if leftover_portion != 0:
            k += 1
            end_coordinate = utr_length
            start_coordinate = end_coordinate - length_limit
            cut_coordinates[k, :] = [k + 1, start_coordinate, end_coordinate]

    else:
        cut_coordinates = np.zeros((1, 3), dtype=int)
        cut_coordinates[0, :] = [1, 0, utr_length]

    return cut_coordinates

def separate_by_rolling_window_v2(fr_name, sequence, length_limit):
    info_string = ''
    cut_intervals = make_intervals_for_rw_separation_v2(len(sequence), length_limit)
    for k in range(cut_intervals.shape[0]):
        curr_list = cut_intervals[k,:]
        info_string += '>%s_part_%d\n' % (fr_name, curr_list[0])
        info_string += sequence[curr_list[1] : curr_list[2]]
        info_string += '\n'
    return info_string

def write_separated(input_dict, output_file, length):
    with open(output_file, 'w') as wf:
        for tr in sorted(list(input_dict.keys())):
            string_to_write = separate_by_rolling_window_v2(tr,
                                                            input_dict[tr],
                                                            length_limit = length)
            wf.write(string_to_write)


def main(raw_args = None):
    args = handler(raw_args)
    fasta_dict = IO.read_fasta(args.f)
    write_separated(fasta_dict, args.o, args.length)


if __name__ == "__main__":
    main()