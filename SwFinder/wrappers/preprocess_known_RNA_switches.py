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





def main(raw_args = None):
    args = handler(raw_args)
    fasta_dict = IO.read_fasta(args.f)
    write_separated(fasta_dict, args.o, args.length)


if __name__ == "__main__":
    main()