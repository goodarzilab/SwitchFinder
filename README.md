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
#### 1. Generate seeds
	Use pyteiser/seeds_generator.py
#### 2. Convert sequences from fasta to binary format
	Use pyteiser/wrappers/binarize_sequences.py
#### 3. Precalculate seed occurence profiles
	Use pyteiser/wrappers/calculate_seed_profiles.py - run on HPC!
#### 4. Preprocess the expression file of interest
	Use either pyteiser/wrappers/preprocess_expression_profile_ensembl.py or pyteiser/wrappers/preprocess_custom_expression_profile.py
#### 5. Calculate MI values for all the seeds
	Use pyteiser/wrappers/calculate_MI_profiles.py - run on HPC!
#### 5a. (optional) Filter possible seed matches with *in silico* RNA folding algorithm
	Use pyteiser/wrappers/filter_profiles_by_folding.py
#### 6. Choose significance thresholds
	Use pyteiser/wrappers/choose_significant_seeds_v3.py - run on HPC!
#### 7. Combine seeds that passed
	Use pyteiser/wrappers/combine_passed_seeds.py
#### 8. Classify seeds by families
	Use pyteiser/wrappers/filter_passed_seeds_with_CMI.py
#### 9. Optimize seeds
	Use pyteiser/wrappers/optimize_seeds_single_chunk.py - run on HPC! You can submit it with pyteiser/wrappers/qsub_optimize_seeds.py
#### 10. Combine optimized seeds
	Use pyteiser/wrappers/combine_optimized_seeds.py

### License
MIT license

### Citing
See the paper

### About SwFinder
SwFinder has been developed in Goodarzi lab at UCSF by Matvei Khoroshkin and Hani Goodarzi