#!/bin/bash

# Build the Docker image
docker build -t eagleshot/switch_finder_image:latest  -t eagleshot/switch_finder_image:1.0 .

## Mount the current directory to the Docker container
# Run the pipeline
docker run -it -v $(pwd):/switchfinder eagleshot/switch_finder_image:latest 
docker run -it -v $(pwd):/code eagleshot/switch_finder_image:latest

## Run the docker image
# docker run -it -v $(pwd):/switchfinder eagleshot/switch_finder_image:latest


## Commands to run the pipeline
# Run the pipeline
python SwitchFinder/wrappers/SwitchFinder_pipeline.py \
        --input_fastafile example_data/example_sequences.fa \
        --out output \
        --temp_folder temp \
        --RNAstructure_path $RNAstructure_path \
        --RNApathfinder_path $RNApathfinder_path

# # Run the new classifier script
python SwitchFinder/wrappers/new_classifier.py \
--input_fastafile example_data/seed_riboswitches.fa \
--temp_folder temp \
--RNAstructure_path $RNAstructure_path \
--RNApathfinder_path $RNApathfinder_path \
--fragment_length 200