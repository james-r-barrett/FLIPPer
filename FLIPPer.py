## check that all modules required to run FLIPPer are installed
import sys
from importlib import util
requirements = ['pandas', 'Bio', 'metapredict', 'cython', 'matplotlib', 'bs4', 'protfasta']

## check that additional packages required for pipeline are installed
## if not detected in installation, ask user if they want to proceed anyway
for requirement in requirements:
    req = util.find_spec(requirement)
    if req is not None:
        continue
    else:
        ignore=(input(requirement+" not detected - proceed anyway? (y/n): "))
        if ignore == 'y':
            print("Ignoring missing packages - errors may be encountered.")
            break
        else:
            print("Exiting. - check package installations and try again")
            sys.exit()

## import modules from python3 standard
import subprocess
import os
import glob
from os import listdir
from os.path import isfile, join

## change import directory scripts folder for import, then import FLIPPer modules
sys.path.insert(1, 'scripts/')
from FLIPPer_lib import *

## define path as current working directory   
PATH = os.getcwd()

## Find all non-directory files in the current directory, then remove FLIPPer default files and OS-specfic hidden files from list
onlyfiles = [f for f in listdir(PATH) if isfile (join(PATH,f))]
package_files = ["FLIPPer.py", "Readme.txt", "Changelog.txt", "desktop.ini", ".DS_Store", ".gitattributes"]
for file in package_files:
	if os.path.isfile(file):
    		onlyfiles.remove(file)

## Define graphics for print outputs
lineenter = "\n" + "=========================================================" + "\n" 

## Check if output folder for files exists
## If it does, remove from files for analysis then continue
for file in onlyfiles:
    if os.path.exists("{}_FLIPPer_outputs".format(file)):
        onlyfiles.remove(file)
        print(lineenter)
        print("Output folder for "+str(file)+" already exists, please rename or remove - skipping.")

## Print blurb
print(lineenter)
print("FLIPPer - Fast Linker Identification Pipeline for Pyrenoids - v1.1 (built 29.02.2024, J. Barrett [james.barrett@york.ac.uk])")
print(lineenter)

## Ask user for inputs to determine how pipeline is run
## Do they want to use custom variables other than default (defined in FLIPPER_lib)
pre_filtering = (input("Would you like to use custom pre-filtering variables? (y/n): "))

## User input section to determine variables for filtering
if pre_filtering == "y":
    print(lineenter + "\n" + "Please set filtering variables..." + "\n")
    pI = float(input("Set pI threshold (Default = 8): "))
    THRatio = float(input("Set Turn/Helix Ratio threshold (Default = 1): "))
    Serine = float(input("Set Serine content threshold (Default = 0.05): "))
    Alanine = float(input("Set Alanine content threshold (Default = 0.01): "))

## Do they want to use custom variables for XSTREAM?
print(lineenter)
xstream_custom = (input("Would you like to use custom XSTREAM variables? (y/n): "))

## User input if so
if xstream_custom == "y":
    print("\nPlease set Xstream variables..." + "\n")
    Copy = input("Set minimum copy number (Default = 3): ")
    Word = input("Set minmum word match (Default = 0.3625): ")
    Consensus = input("Set consensus match (Default = 0.4): ")
    Gaps = input("Set maximum gaps in repeats (Default = 50): ")
    minPeriod = input("Set minimum repeat period (Default = 20): ")
    maxPeriod = input("Set maximum repeat period (Default = 120): ")
    Coverage = input("Set sequence proportion repeats cover (0.4 for fusions, 0.75 for full linkers): ")

## Do they want to filter sequences after XSTREAM to only include sequences with repeats that contain aromatic/electrostatic residues?
print(lineenter)
xstream_positive = (input("Would you like to use custom filtering paramaters for sequences after XSTREAM? (Default = 1 aromatic, 2 electrostatic) (y/n): "))
## If yes to above, specify number of electrosatic and aromatic resiudes
if xstream_positive == "y":
    print("\n" + "Please set post-XSTREAM filtering variables..." + "\n")
    Aromatic = float(input("Set minimum number of aromatic resiudes in repeat region (Default = 1): "))
    Electrostatic = float(input("Set minmum number of electrostatic residues in repeat region (Default = 2): "))

## Do they want to keep the sequence analysis of the full input file after analysis and filtering, otherwise keep the candidates
print(lineenter)
full_output = (input("Would you like to keep sequence analysis of all input sequences? - large file, n recommended (y/n): "))

## Do they want to filter after metapredict analysis?
print(lineenter)
metapredict_filter = (input("Would you like to use a custom disorder filter after metapredict analysis? (y/n): "))
if metapredict_filter == 'y':
        metapredict_filter_value = float(input("Threshold for IUPred filtering? (Default = 50): "))
else:
        metapredict_filter_value = 50

## Do they want to output metapredict plots for each of the XSTREAM sequences?
print(lineenter)
metapredict_plot = (input("Would you like to plot metapredict/pLDDT profiles for candidate sequences? - y recommended (y/n): "))

## for each file in onlyfiles without output directory already existing
for file in onlyfiles:
    print(file)
    ## validate if the file is fasta, if true, complete pipeline
    while validate_fasta(file):   
        ## make directories for outputs
        destination_folder="{}_FLIPPer_outputs".format(file)
        directory = "{}_FLIPPer_outputs/metapredict_plots".format(file)
        os.mkdir(destination_folder)
        os.makedirs(directory)
        
        ## run analysis and filtering module
        analysis_and_filtering(file, pI, THRatio, Serine, Alanine, Aromatic, Electrostatic, full_output)
        
        ## run XSTREAM using filtered sequences
        subprocess.call(["java", "-Xmx1000m", "-Xms1000m", "-jar", "scripts/xstream.jar", "CandidateSequences_Temp.FASTA", "-e"+str(Copy), "-i"+str(Word), "-I" +str(Consensus), "-g" +str(Gaps), "-m" +str(minPeriod), "-x" +str(maxPeriod), "-a"+str(file), "-t1", "-T" +str(Coverage)])
        
        ## find html file which is output from XSTREAM, extract sequences, then run metapredict module against
        for f in os.listdir(PATH):
            if f.endswith("_2.html"):
                ## extract sequence IDs, then retrieve sequences from input file and output temp file for iupred/filtering
                ## if XSTREAM filtering used, repeat region extract from _2.html file, filtered with input variables
                ## then sequence IDs used to extract sequences of filtered sequences from input file
                xstream_extract2(f, file, Aromatic, Electrostatic)
                ## remove XSTREAM files before second run
                for x in glob.glob("XSTREAM*"):
                    os.remove(x)
                ## run XSTREAM again on filtered sequences
                ## not ideal in terms of speed, but cleans up the outputs
                subprocess.call(['java', '-Xmx1000m', '-Xms1000m', '-jar', 'scripts/xstream.jar', 'Temp_xstream_filtered.fasta', "-e"+str(Copy), "-i"+str(Word), "-I" +str(Consensus), "-g" +str(Gaps), "-m" +str(minPeriod), "-x" +str(maxPeriod), "-a"+str(file), "-t1", "-T" +str(Coverage)], stdout=subprocess.DEVNULL)   
                print(lineenter)
                print("XSTREAM analysis of filtered proteins finished!")
                print(lineenter)

                ## run metapredict module against filtered sequences
                metapredict_htp('Temp_xstream_filtered.fasta',directory, metapredict_plot, metapredict_filter_value)

                for x in glob.glob("XSTREAM*"):
                    os.remove(x)
                ## run XSTREAM again on filtered sequences
                ## not ideal in terms of speed, but cleans up the outputs
                subprocess.call(['java', '-Xmx1000m', '-Xms1000m', '-jar', 'scripts/xstream.jar', 'candidate_sequences.fasta', "-e"+str(Copy), "-i"+str(Word), "-I" +str(Consensus), "-g" +str(Gaps), "-m" +str(minPeriod), "-x" +str(maxPeriod), "-a"+str(file), "-t1", "-T" +str(Coverage)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) 

        print(lineenter)
        
        ## pass variables from input to output_variables module to write variables file
        output_variables(file, pI, THRatio, Serine, Alanine, Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage, Type, Aromatic, Electrostatic, metapredict_filter_value)
        
        ## cleanup files into _outputs folder
        print("Cleaning up files from "+str(file)+ " analysis.")
        os.remove("CandidateSequences_Temp.FASTA")
        os.remove("Temp_xstream_filtered.fasta")
        os.remove("Temp_Xstream_positives.fasta")
        for z in glob.glob("XSTREAM*"):
            move_files(z, destination_folder, PATH)      
        move_files("sequence_analysis", destination_folder, PATH)
        move_files(file, destination_folder, PATH)
        move_files(file+"_variables.txt", destination_folder, PATH)
        move_files("candidate_sequences.fasta", destination_folder, PATH)
        move_files("candidate_sequences.csv", destination_folder, PATH)
        
        ## end of analysis run
        print("Analysis of "+str(file)+" finished!") 
        print(lineenter)
        
        ## Update onlyfiles list so loop doesn't search for file just analysed
        onlyfiles = [f for f in listdir(PATH) if isfile (join(PATH,f))]
        break
    
    ## if file is not in fasta format, skip it and return this
    else:
        print(lineenter)
        print(file + " is not FASTA format - skipping.")
        print(lineenter) 

print("ALL FINISHED!")

















