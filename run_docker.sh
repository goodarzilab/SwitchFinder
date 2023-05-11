#!/bin/bash

# Build the Docker image
docker build -t switch_finder_image .

## Mount the current directory to the Docker container
# Run the pipeline
docker -v {pwd}:/SwitchFinder run switch_finder_image python \
        /SwitchFinder/SwitchFinder/wrappers/SwitchFinder_pipeline.py \
        --input_fastafile /SwitchFinder/example_data/example_sequences.fa \
        --out /khayyam/khorms/temp/SwitchFinder/out \
        --temp_folder /khayyam/khorms/temp/SwitchFinder/temp \
        --RNAstructure_path $RNAstructure_path \
        --RNApathfinder_path $RNApathfinder_path

# Run the new classifier script
docker -v {pwd}:/SwitchFinder run switch_finder_image \
python /SwitchFinder/SwitchFinder/wrappers/new_classifier.py \
--input_fastafile /khayyam/khorms/temp/SwitchFinder/seed_riboswitches.fa \
--temp_folder /khayyam/khorms/temp/SwitchFinder/temp \
--RNAstructure_path $RNAstructure_path \
--RNApathfinder_path $RNApathfinder_path \
--fragment_length 200
