import os
import argparse

import sys

sys.path.append('/avicenna/khorms/programs/SwFinder/SwFinder')
sys.path.append('/avicenna/khorms/programs/SwFinder/SwFinder/wrappers')

import chop_sequences as chop_sequences

import IO as IO

def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fastafile", help="the fasta file with the target sequences", type=str)
    parser.add_argument("-out", help="output folder", type=str)
    parser.add_argument("--temp_folder", help="folder for temporary files", type=str)
    parser.add_argument("--fragment_length", help="fragment length", type=int)


    parser.set_defaults(
                        input_fastafile = "/avicenna/khorms/programs/SwFinder/example_data/example_sequences.fa",
                        out = "/avicenna/khorms/temp/SwFinder/out",
                        temp_folder = "/avicenna/khorms/temp/SwFinder/temp",
                        fragment_length = 186,
                        )
    args = parser.parse_args()
    return args


def generate_filenames(args):
    # list folders
    temp_folder = args.temp_folder
    out_folder = args.out
    inputs_folder = os.path.join(temp_folder, "inputs")
    interm_folder = os.path.join(temp_folder, "interm")

    # list filenames
    chopped_sequences_filename = os.path.join(interm_folder, "chopped_sequences.fa")

    # make directories
    IO.create_folder(temp_folder)
    IO.create_folder(inputs_folder)
    IO.create_folder(interm_folder)
    IO.create_folder(out_folder)

    filenames_dict = {
        "chopped_sequences" : chopped_sequences_filename,
    }

    return filenames_dict


def main():
    args = handler()
    filenames_dict = generate_filenames(args)

    chop_sequences_args = ["-f", args.input_fastafile,
                           "-o", filenames_dict['chopped_sequences'],
                           "--length", str(args.fragment_length)]
    chop_sequences.main(chop_sequences_args)


if __name__ == '__main__':
    main()