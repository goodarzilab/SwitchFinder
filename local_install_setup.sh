#!/bin/bash

# Make new directory if it doesn't exist
mkdir -p programs

# Copy files from remote server
scp -r user@avicenna:/khorms/programs/RNApathfinder ./programs
scp -r user@avicenna:/khorms/programs/RNAstructure/ ./programs