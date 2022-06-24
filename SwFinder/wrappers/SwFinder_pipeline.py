import os
import argparse
import shutil

import sys

sys.path.append('/avicenna/khorms/programs/SwFinder/SwFinder')
sys.path.append('/avicenna/khorms/programs/SwFinder/SwFinder/wrappers')

import chop_sequences
import find_mutually_exclusive_stems
import fold_mutually_exclusive_structures
import calculate_energy_barriers
import apply_classifier
import generate_mutations

import IO

def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_fastafile", help="the fasta file with the target sequences", type=str)
    parser.add_argument("-out", help="output folder", type=str)
    parser.add_argument("--temp_folder", help="folder for temporary files", type=str)
    parser.add_argument("--fragment_length", help="fragment length", type=int)
    parser.add_argument("--RNAstructure_path", help="path to RNAstructure program", type=str)
    parser.add_argument("--RNApathfinder_path", help="path to RNApathfinder program", type=str)
    parser.add_argument("--num_processes", help="number of processes to run in parralel", type=int)
    parser.add_argument("--n_iterations", help="number of iterations for the mutation generator script", type=int)
    parser.add_argument("--n_mutations", help="number of mutations to generate for each conformation", type=int)


    parser.set_defaults(
                        input_fastafile = "/avicenna/khorms/programs/SwFinder/example_data/example_sequences.fa",
                        out = "/avicenna/khorms/temp/SwFinder/out",
                        temp_folder = "/avicenna/khorms/temp/SwFinder/temp",
                        RNAstructure_path='/avicenna/khorms/programs/RNAstructure',
                        RNApathfinder_path='/avicenna/khorms/programs/RNApathfinder',
                        num_processes = 20,
                        fragment_length = 186,
                        n_iterations = 100,
                        n_mutations = 1,
                        )
    args = parser.parse_args()
    return args


def generate_filenames(args):
    # list folders
    temp_folder = args.temp_folder
    out_folder = args.out
    inputs_folder = os.path.join(temp_folder, "inputs")
    interm_folder = os.path.join(temp_folder, "interm")
    MIBP_temp_folder = os.path.join(interm_folder, "MIBP")

    # list filenames
    chopped_sequences_filename = os.path.join(interm_folder, "chopped_sequences.fa")
    mutually_exclusive_stems_folder = os.path.join(interm_folder, "mutually_exclusive_stems_folder")
    mutually_exclusive_stems_filename = os.path.join(mutually_exclusive_stems_folder, 'output.txt')
    mutually_exclusive_conformations_filename = os.path.join(interm_folder, 'mutually_exclusive_conformations.txt')
    energy_barriers_filename = os.path.join(interm_folder, 'energy_barriers.txt')
    classifier_output_filename = os.path.join(out_folder, "RNA_switch_scores.txt")
    pipeline_text_output_filename = os.path.join(out_folder, "RNA_switch_structures_full.txt")
    pipeline_text_output_short_filename = os.path.join(out_folder, "RNA_switch_structures.txt")
    generate_mutations_folder = os.path.join(interm_folder, "generate_mutations_folder")
    generate_mutations_filename = os.path.join(generate_mutations_folder, 'perturbations.txt')
    mutations_final_filename = os.path.join(out_folder, "generated_mutations.txt")


    # make directories
    IO.create_folder(temp_folder)
    IO.create_folder(inputs_folder)
    IO.create_folder(interm_folder)
    IO.create_folder(out_folder)
    IO.create_folder(MIBP_temp_folder)
    IO.create_folder(mutually_exclusive_stems_folder)
    IO.create_folder(generate_mutations_folder)

    filenames_dict = {
        "temp_folder" : temp_folder,
        "MIBP_temp_folder" : MIBP_temp_folder,
        "chopped_sequences" : chopped_sequences_filename,
        "mutually_exclusive_stems_folder" : mutually_exclusive_stems_folder,
        "mutually_exclusive_stems_filename" : mutually_exclusive_stems_filename,
        "mutually_exclusive_conformations_filename" : mutually_exclusive_conformations_filename,
        "energy_barriers_filename" : energy_barriers_filename,
        "classifier_output_filename" : classifier_output_filename,
        "pipeline_text_output_filename" : pipeline_text_output_filename,
        "pipeline_text_output_short_filename" : pipeline_text_output_short_filename,
        "generate_mutations_folder" : generate_mutations_folder,
        "generate_mutations_filename" : generate_mutations_filename,
        "mutations_final_filename" : mutations_final_filename
    }

    return filenames_dict

def clean_up_afterwards(filenames_dict):
    shutil.copy(filenames_dict["generate_mutations_filename"], filenames_dict["mutations_final_filename"])
    shutil.rmtree(filenames_dict["temp_folder"], ignore_errors=True)
    os.remove(filenames_dict["pipeline_text_output_filename"])
    IO.create_folder(filenames_dict["temp_folder"])


def main():
    args = handler()
    filenames_dict = generate_filenames(args)

    chop_sequences_args = [
        "-f", args.input_fastafile,
        "-o", filenames_dict['chopped_sequences'],
        "--length", str(args.fragment_length)
    ]
    find_mutually_exclusive_stems_args = [
        "-f", filenames_dict['chopped_sequences'],
        "--temp_files_folder", filenames_dict['MIBP_temp_folder'],
        "-o", filenames_dict['mutually_exclusive_stems_folder'],
        "--RNAstructure_path", args.RNAstructure_path,
        "--num_processes", str(args.num_processes),
    ]
    fold_mutually_exclusive_structures_args = \
        ["-f", filenames_dict['chopped_sequences'],
        "--temp_files_folder", filenames_dict['MIBP_temp_folder'],
        "-s", filenames_dict['mutually_exclusive_stems_filename'],
        "-o", filenames_dict['mutually_exclusive_conformations_filename'],
        "--RNAstructure_path", args.RNAstructure_path,
        "--num_processes", str(args.num_processes),
    ]
    calculate_energy_barriers_args = [
        "--dotbracket", filenames_dict['mutually_exclusive_conformations_filename'],
        "-o", filenames_dict['energy_barriers_filename'],
        "--temp_files_folder", filenames_dict['MIBP_temp_folder'],
        "--num_processes", str(args.num_processes),
        "--path_rnapathfinder", args.RNApathfinder_path,
    ]
    apply_classifier_args = [
        "--energies_filename", filenames_dict['energy_barriers_filename'],
        "--dataframe_output", filenames_dict['classifier_output_filename'],
        "--text_output", filenames_dict['pipeline_text_output_filename'],
        "--text_output_short", filenames_dict["pipeline_text_output_short_filename"]
    ]
    generate_mutations_args = [
        "-f", filenames_dict['chopped_sequences'],
        "-s", filenames_dict['mutually_exclusive_stems_filename'],
        "-o", filenames_dict["generate_mutations_folder"],
        "--n_iterations", str(args.n_iterations),
        "--n_loops_desired", str(args.n_mutations),
    ]

    chop_sequences.main(chop_sequences_args)
    find_mutually_exclusive_stems.main(find_mutually_exclusive_stems_args)
    fold_mutually_exclusive_structures.main(fold_mutually_exclusive_structures_args)
    calculate_energy_barriers.main(calculate_energy_barriers_args)
    apply_classifier.main(apply_classifier_args)
    generate_mutations.main(generate_mutations_args)
    clean_up_afterwards(filenames_dict)



if __name__ == '__main__':
    main()