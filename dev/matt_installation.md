# create environment
mamba create -n switch_finder_env python 
mamba install -y -c conda-forge pandas numpy scikit-learn statsmodels
mamba install -y -c bioconda viennarna
# copy over pre-compiled RNApathfinder and RNAstructure from avicenna
scp -r /avicenna/khorms/programs/RNApathfinder khayyam:/khayyam/khorms/programs
scp -r /avicenna/khorms/programs/RNAstructure/ khayyam:/khayyam/khorms/programs
# run the pipeline
python /khayyam/khorms/programs/SwitchFinder/SwitchFinder/wrappers/SwitchFinder_pipeline.py --input_fastafile /khayyam/khorms/programs/SwitchFinder/example_data/example_sequences.fa --out /khayyam/khorms/temp/SwitchFinder/out --temp_folder /khayyam/khorms/temp/SwitchFinder/temp --RNAstructure_path /khayyam/khorms/programs/RNAstructure --RNApathfinder_path /khayyam/khorms/programs/RNApathfinder
# make sure the "new classifier" script also works
python /khayyam/khorms/programs/SwitchFinder/SwitchFinder/wrappers/new_classifier.py --input_fastafile /khayyam/khorms/temp/SwitchFinder/seed_riboswitches.fa --temp_folder /khayyam/khorms/temp/SwitchFinder/temp --RNAstructure_path /khayyam/khorms/programs/RNAstructure --RNApathfinder_path /khayyam/khorms/programs/RNApathfinder --fragment_length 200