#!/bin/bash

# Build the Docker image
docker build -t switch_finder_image .

## Mount the current directory to the Docker container
# Run the pipeline
docker -v {pwd}:/SwitchFinder run switch_finder_image 

## Run the docker image
# docker run -it -v $(pwd):/switchfinder switch_finder_image


## Commands to run the pipeline
# Run the pipeline
python SwitchFinder/wrappers/SwitchFinder_pipeline.py \
        --input_fastafile example_data/example_sequences.fa \
        --out output \
        --temp_folder temp \
        --RNAstructure_path $RNAstructure_path \
        --RNApathfinder_path $RNApathfinder_path

# # Run the new classifier script
# docker -v {pwd}:/SwitchFinder run switch_finder_image python SwitchFinder/wrappers/new_classifier.py \
# --input_fastafile example_data/seed_riboswitches.fa \
# --temp_folder temp \
# --RNAstructure_path $RNAstructure_path \
# --RNApathfinder_path $RNApathfinder_path \
# --fragment_length 200
