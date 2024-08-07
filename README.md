# SwitchFinder
SwitchFinder is a computational tool designed for the systematic discovery of RNA structural switches within transcriptomes. It operates by analyzing RNA sequences to predict potential RNA switches and their two mutually exclusive folding conformations. Utilizing a biophysical approach based on a Boltzmann equilibrium probability distribution, SwitchFinder identifies RNA switches in a family-agnostic manner, prioritizing sequences with features indicative of RNA switches. See [Khoroshkin et al, 2024](https://www.nature.com/articles/s41592-024-02335-1)

## Installation
### Using Docker (recommended)
Prerequisites:
- Working installation of Docker. Follow the official installation guide.
- Up to date installation of pip. Pip can be upgraded with the command pip install --upgrade pip setuptools wheel
Steps:
1. Pull a docker container image with the following command (Running the command might require “sudo” depending on your user privileges)
	```
	docker pull goodarzilaborder/switch_finder_image:latest
	```
2. Launch the Docker image
	```
	docker run -it -v $(pwd):/switchfinder goodarzilaborder/switch_finder_image:latest
	```
3. Install the Python code inside the Docker container
	```
	pip install git+https://github.com/goodarzilab/SwitchFinder.git
	```
4. Test the installation

    1. Download example data (within the Docker container):
	```
	wget https://raw.githubusercontent.com/goodarzilab/SwitchFinder/main/example_data/example_sequences.fa
	```
    2. Run the pipeline (within the Docker container)
	```
	SwitchFinder_pipeline \
	        --input_fastafile example_sequences.fa \
	        --out output \
	        --temp_folder temp \
	        --RNAstructure_path $RNAstructure_path \
	        --RNApathfinder_path $RNApathfinder_path
	```
    3. Move the output files to the desired location (outside of Docker container)



### Alternative installation (individual packages)
SwitchFinder is meant to be run on Linux. The software has not been tested on Windows or Mac OS
1. Install RNApathfinder
	- Download the RNApathfinder software from this repo to a specified folder ([RNAPathFinder.tar](https://github.com/goodarzilab/SwitchFinder/blob/main/RNAPathFinder.tar)). Please cite I. Dotú, W.A. Lorenz, P. Van Hentenryck, P. Clote. Nucleic Acids Res. 2010 Mar 1;38(5):1711-22.
	- unpack the archive with the command `tar -xf RNAPathFinder.tar`
 	- compile as a standard C software
	- pass the address of this folder (titled srcTABU) as an argument to SwitchFinder scripts
2. Install RNAstructure
	- Download the RNAstructure software from [here](http://rna.urmc.rochester.edu/Releases/current/RNAstructureForLinux.tgz) or [here](https://rna.urmc.rochester.edu/RNAstructure.html)
	- unpack the archive with the command `tar -zxf RNAstructureForLinux.tgz`
	- pass the address of this folder as an argument to SwitchFinder scripts	
3. Install SwitchFinder:
```
pip install git+https://github.com/goodarzilab/SwitchFinder.git
```

### Usage of the automatic pipeline
First, download the example input file with [sequences](https://github.com/goodarzilab/SwitchFinder/blob/main/example_data/example_sequences.fa)<br>
Next, make sure that the Python requirements, along with RNApathfinder and RNAstructure, are installed
Then, launch the pipeline with the command: <br>
```
python SwitchFinder_pipeline.py --input_fastafile <path to example_sequences.fa> -out <path to the output folder> --temp_folder <path to the folder for temporary files> --RNAstructure_path <path to the RNAstructure installation directory> --RNApathfinder_path <path to the RNApathfinder installation directory>
```

Input files:
- Necessary:
	- `input_fastafile`: a fasta file with RNA sequences of interest
- Optional:
	- user can include a file with RNA structure probing data (SHAPE or DMS-seq) to guide the possible match selection. There is no commonly used standard format for SHAPE RNA reactivity data. SHAPE file provided by user should be a pickled dictionary, where the names are the names of individual fragments, and the values are single-dimensional numpy arrays containing the normalized SHAPE/DMSseq values. The positions where the data is absent are masked by negative values (for example, -999). SHAPE file can be provided to the `find_mutually_exclusive_stems.py` and `fold_mutually_exclusive_structures.py` scripts with the `--shape_profile` argument

Output files:
The pipeline generates three files:
- `RNA_switch_scores.txt`: a table with scores for predicted RNA switches. Each row corresponds to a single predicted switch. The two columns contain the names and the scores (higher score corresponds to higher probability of RNA switch)
- `RNA_switch_structures.txt`: a text file with predicted RNA switches and the mutually exclusive conformations. The descriptions of individual switches are separated by `$$$`. The switches are listed in the same order as in `RNA_switch_scores.txt`; even low scoring switches are listed. We leave it up to the user to choose the top N predictions. The description of an individual switch includes the name, the sequence, the secondary structures of the two mutually exclusive conformations, and the substructures that are common for the two conformations
- `generated_mutations.txt`: a text file with predicted sequence mutations that shift the equilibrium between the two mutually exclusive confromations. The descriptions of individual switches are separated by `$$$`. All of the predicted switches (even the low scoring ones) are listed. The description of an individual switch includes: the name, the sequence, the coordinates of the key stem loops of the two mutually exclusive confromations, the predicted sequence mutations of four different categories. The categories are: (1) major_strengthen: strengthen base pairing of the conformation 1; second_strengthen: strengthen base pairing of the conformation 2; second_weaken: disrupt the base pairing of the conformation 2; major_weaken: disrupt the base pairing of the conformation 1 


### Usage of pipeline in manual step-by-step mode
#### 1. Split the sequence to fragments of the same length
	Use chop_sequences

#### 2. Identify the local minima of the RNA folding landscape
	Use find_mutually_exclusive_stems

#### 3. Fold the mutually exclusive structures
	Use fold_mutually_exclusive_structures

#### 4. Calculate the energies required for transition between the two conformations
	Use calculate_energy_barriers

#### 5. Predict which fragments are likely to be RNA switches
	Use apply_classifier

#### 6. Generate sequence mutations that shift the equilibrium between the two mutually exclusive confromations
	Use cgenerate_mutations

### Generating your own classifier
At the step 5, we assign scores to the individual predicted RNA switches. We use scores from a classifier that was pre-trained on a set of known bacterial riboswitches, downloaded from [Rfam](https://rfam.xfam.org/). If you wish to pre-train your own classifier, you may use a script named `new_classifier.py`. As the input file, please provide a fasta file containing only the known riboswitches. As an example, we provide a fasta [file](https://github.com/goodarzilab/SwitchFinder/blob/main/example_data/seed_riboswitches.fa) with the sequences of known riboswitches downloaded from [Rfam](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT). You should specify the parameter `--fragment_length`; all the sequences longer than the value you specify will be ignored. We recommend the value of 200, since most known RNA switches are between 0-200 nt long. The `new_classifier.py` script will output 3 parameters for a classifier; to apply them to the set of sequences of interest, pass them to the `SwitchFinder_pipeline.py` pipeline as parameters `--loop_energies_coefficient`, `--barrier_heights_coefficient`, `--intercept`. The parameters you have to specify:
```
new_classifier --input_fastafile <path to known RNA switch sequences.fa> --temp_folder <path to the folder for temporary files> --RNAstructure_path <path to the RNAstructure installation directory> --RNApathfinder_path <path to the RNApathfinder installation directory> --fragment_length <desired length limit>
```

### License
MIT license

### Citing
[Khoroshkin et al, 2024](https://www.nature.com/articles/s41592-024-02335-1)

### About SwitchFinder
SwitchFinder has been developed in Goodarzi lab at UCSF by Matvei Khoroshkin and Hani Goodarzi
