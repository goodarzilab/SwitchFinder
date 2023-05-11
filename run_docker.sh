#!/bin/bash

# Build the Docker image
docker build -t switch_finder_image .

# Run the pipeline
docker run switch_finder_image python /khayyam/khorms/programs/SwitchFinder/SwitchFinder/wrappers/SwitchFinder_pipeline.py --input_fastafile /khayyam/khorms/programs/SwitchFinder/example_data/example_sequences.fa --out /khayyam/khorms/temp/SwitchFinder/out --temp_folder /khayyam/khorms/temp/SwitchFinder/temp --RNAstructure_path /khayyam/khorms/programs/RNAstructure --RNApathfinder_path /khayyam/khorms/programs/RNApathfinder

# Run the new classifier script
docker run switch_finder_image python /khayyam/khorms/programs/SwitchFinder/SwitchFinder/wrappers/new_classifier.py --input_fastafile /khayyam/khorms/temp/SwitchFinder/seed_riboswitches.fa --temp_folder /khayyam/khorms/temp/SwitchFinder/temp --RNAstructure_path /khayyam/khorms/programs/RNAstructure --RNApathfinder_path /khayyam/khorms/programs/RNApathfinder --fragment_length 200
