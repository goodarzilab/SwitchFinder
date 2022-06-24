# SwFinder
Software for predicting RNA switches

### Introduction
SwFinder 

### Installation
Install with pip: <br>
First, install the requirements:
```
pip install numba==0.50.1 numpy>=1.19.1 pandas>=1.1.1
```
Then, install pyteiser:
```
pip install pyteiser
```
### Usage of the automatic pipeline
First, download the three example input files: [sequences](https://github.com/goodarzilab/pyteiser/raw/master/example_data/test_seqs.fa), [measurements](https://github.com/goodarzilab/pyteiser/raw/master/example_data/test_expression_values.csv) and [seeds](https://github.com/goodarzilab/pyteiser/raw/master/example_data/test_seeds_20.bin) <br>
Then, create folders for temporary files and for output files (you can do it with `mkdir temp` and `mkdir out`)<br>
Finally, launch the pipeline with the command: <br>
```
pyteiser_pipeline --rna_fastafile <path to test_seqs.fa> --exp_values_file <path to test_expression_values.csv> --seeds_file <path to test_seeds_20.bin> --temp_folder <path to the folder for temporary files> --out <path to the output folder>
```


### Usage of pipeline in manual step-by-step mode
#### 1. Split the sequence to fragments of the same length
chop_sequences.py

#### 2. Identify the conflicting local minima
find_mutually_exclusive_stems.py

#### 2. Fold the mutually exclusive structures
	Use SwFinder/MIBP_fragments_stem_finder_SHAPE_parralel.py

#### 2. Calculate the energies required for transition between the two conformations
	Use SwFinder/energy_differences_v7_RNApathfinder.py

#### 2. Train a classifier
	Use SwFinder/MIBP_fragments_stem_finder_SHAPE_parralel.py	

#### 2. Apply a classifier
	Use SwFinder/convert_MIBP_output_for_logit_classifier.ipynb -> convert_energy_calculations_to_df_v3
	some simple Python script?	



### License
MIT license

### Citing
See the paper

### About SwFinder
SwFinder has been developed in Goodarzi lab at UCSF by Matvei Khoroshkin and Hani Goodarzi