# SwFinder
Software for predicting RNA switches at a genome-wide scale

### Introduction
SwFinder explicityly models 

### Installation
SwFinder is meant to be run on Linux. The software has not been tested on Windows or Mac OS
First, install the required Python packages:
```
pip install numpy pandas sklearn
```
Then, install the required external software:
- RNApathfinder
	- Download the RNApathfinder software from [here](http://bioinformatics.bc.edu/clotelab/RNApathfinder/srcTABU.tgz) to a specified folder
	- unpack the archive with the command `tar -zxf srcTABU.tgz`
	- pass the address of this folder as an argument to SwFinder scripts
- RNAstructure
	- Download the RNAstructure software from [here](http://rna.urmc.rochester.edu/Releases/current/RNAstructureForLinux.tgz) or [here](https://rna.urmc.rochester.edu/RNAstructure.html)
	- unpack the archive with the command `tar -zxf RNAstructureForLinux.tgz`
	- pass the address of this folder as an argument to SwFinder scripts	

Then, install SwFinder:
```
pip install SwFinder
```
### Usage of the automatic pipeline
First, download the example input file with [sequences](https://github.com/goodarzilab/SwFinder/blob/main/example_data/example_sequences.fa)<br>
Next, make sure that the Python requirements, along with RNApathfinder and RNAstructure, are installed
Then, launch the pipeline with the command: <br>
```
python SwFinder_pipeline.py --input_fastafile <path to example_sequences.fa> -out <path to the output folder> --temp_folder <path to the folder for temporary files> --RNAstructure_path <path to the RNAstructure installation directory> --RNApathfinder_path <path to the RNApathfinder installation directory>
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
	Use chop_sequences.py

#### 2. Identify the local minima of the RNA folding landscape
	Use find_mutually_exclusive_stems.py

#### 3. Fold the mutually exclusive structures
	Use fold_mutually_exclusive_structures.py

#### 4. Calculate the energies required for transition between the two conformations
	Use calculate_energy_barriers.py

#### 5. Predict which fragments are likely to be RNA switches
	Use apply_classifier.py

#### 6. Generate sequence mutations that shift the equilibrium between the two mutually exclusive confromations
	Use cgenerate_mutations.py

### Generating your own classifier
At the step 5, we assign scores to the individual predicted RNA switches. We use scores from a classifier that was pre-trained on a set of known bacterial riboswitches, downloaded from [Rfam](https://rfam.xfam.org/). If you wish to pre-train your own classifier, you may use a script named `new_classifier.py`. As the input file, please provide a fasta file containing only the known riboswitches. You should specify the parameter `--fragment_length`; all the sequences longer than the value you specify will be ignored. We recommend the value of 200, since most known RNA switches are between 0-200 nt long. The `new_classifier.py` script will output 3 parameters for a classifier; to apply them to the set of sequences of interest, pass them to the `SwFinder_pipeline.py` pipeline as parameters `--loop_energies_coefficient`, `--barrier_heights_coefficient`, `--intercept`. The parameters you have to specify:
```
python new_classifier.py --input_fastafile <path to known RNA switch sequences.fa> -out <path to the output folder> --temp_folder <path to the folder for temporary files> --RNAstructure_path <path to the RNAstructure installation directory> --RNApathfinder_path <path to the RNApathfinder installation directory> --fragment_length <desired length limit>
```

### License
MIT license

### Citing
See the paper

### About SwFinder
SwFinder has been developed in Goodarzi lab at UCSF by Matvei Khoroshkin and Hani Goodarzi