## check that all modules required to run FLIPPer are installed
import sys
import os
import argparse
from importlib import util

## resolve paths relative to this script's own location, not the caller's working directory -
## this lets FLIPPer be installed once and run against data in any directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

## define command line interface, used both for --help and for optional non-interactive mode
def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="FLIPPer - Fast Linker Identification Pipeline for Pyrenoids"
    )
    parser.add_argument("--engine", choices=["xstream", "detectrepeats"], default=None, help="Repeat-detection engine to use (default: xstream when omitted; prompted for under --interactive)")
    parser.add_argument("--input-dir", default=".", help="Directory containing the FASTA files to process (default: current directory)")
    parser.add_argument("--characterize", metavar="FASTA", help="Characterize a FASTA file of known/reference target sequences (e.g. known pyrenoid linkers) - reports their pI, repeat structure and disorder, and suggests FLIPPer search parameters from them. Runs standalone and exits without processing --input-dir.")
    parser.add_argument("--refilter", metavar="OUTPUT_DIR", help="[--engine detectrepeats] Re-apply --min-period/--max-period/--min-copies/--coverage/--aromatic/--electrostatic/--metapredict-filter-value to an existing <file>_FLIPPer_outputs folder from a previous detectrepeats run, without re-running detection on the full input, and rebuild the candidate report. Writes into OUTPUT_DIR/refiltered/. --min-score cannot be changed this way (it's a detection-time cutoff, not a post-hoc filter) - re-run FLIPPer fully to change it. Runs standalone and exits without processing --input-dir.")
    parser.add_argument("--interactive", action="store_true", help="Prompt for each parameter instead of using the flags below/their documented defaults - press Enter at any prompt to keep the default shown in brackets")
    parser.add_argument("--non-interactive", action="store_true", help=argparse.SUPPRESS)  ## deprecated no-op: running with no flags now behaves this way by default
    parser.add_argument("--pi", type=float, default=8.0, help="pI threshold (default: 8)")
    parser.add_argument("--pi-direction", choices=["min", "max"], default="min", help="Whether --pi keeps proteins with pI >= threshold ('min', default - for basic/Arg-Lys-rich linkers like EPYC1/CsLinker) or pI <= threshold ('max' - for acidic target families)")
    parser.add_argument("--th-ratio", type=float, default=1.0, help="Turn/Helix ratio threshold (default: 1)")
    parser.add_argument("--serine", type=float, default=0.05, help="Serine content threshold (default: 0.05)")
    parser.add_argument("--alanine", type=float, default=0.01, help="Alanine content threshold (default: 0.01)")
    parser.add_argument("--copy", default="3", help="[--engine xstream] XSTREAM minimum copy number (default: 3)")
    parser.add_argument("--word", default="0.3625", help="[--engine xstream] XSTREAM minimum word match (default: 0.3625)")
    parser.add_argument("--consensus", default="0.4", help="[--engine xstream] XSTREAM consensus match (default: 0.4)")
    parser.add_argument("--gaps", default="55", help="[--engine xstream] XSTREAM maximum gaps in repeats (default: 55)")
    parser.add_argument("--min-score", default="8", help="[--engine detectrepeats] DetectRepeats minimum repeat significance score (default: 8) - a detection-time cutoff, not changeable via --refilter")
    parser.add_argument("--min-copies", default="3", help="[--engine detectrepeats] Minimum tandem repeat copy number (default: 3)")
    parser.add_argument("--min-period", default="20", help="Minimum repeat period/unit length in aa (default: 20)")
    parser.add_argument("--max-period", default="120", help="Maximum repeat period/unit length in aa (default: 120)")
    parser.add_argument("--coverage", default="0.75", help="Minimum fraction of the sequence the repeat region must cover - 0.4 for fusions, 0.75 for full linkers (default: 0.75)")
    parser.add_argument("--aromatic", type=float, default=1.0, help="Minimum aromatic residues in repeat region (default: 1)")
    parser.add_argument("--electrostatic", type=float, default=2.0, help="Minimum electrostatic residues in repeat region (default: 2)")
    parser.add_argument("--full-output", action="store_true", help="Keep the full sequence analysis of all input sequences (large file)")
    parser.add_argument("--metapredict-filter-value", type=float, default=50.0, help="Minimum metapredict disorder percentage required within the identified repeat region (default: 50)")
    parser.add_argument("--plots", dest="metapredict_plot", action="store_true", default=True, help="Plot metapredict/pLDDT profiles for candidate sequences (default: on)")
    parser.add_argument("--no-plots", dest="metapredict_plot", action="store_false", help="Disable metapredict/pLDDT plotting")
    return parser

args = build_arg_parser().parse_args()

## module to check a list of required packages is importable, prompting to proceed anyway (or
## exiting outright in --non-interactive mode) if one isn't - used both for the packages every
## run needs and, once the engine is known, for that engine's own extra requirements
def check_requirements(requirements, non_interactive):
    for requirement in requirements:
        req = util.find_spec(requirement)
        if req is not None:
            continue
        if non_interactive:
            print(requirement + " not detected. Exiting.")
            sys.exit(1)
        ignore = (input(requirement+" not detected - proceed anyway? (y/n): "))
        if ignore == 'y':
            print("Ignoring missing packages - errors may be encountered.")
            break
        else:
            print("Exiting. - check package installations and try again")
            sys.exit()

## check that packages required by every run (regardless of engine) are installed - engine-specific
## extras (e.g. bs4 for --engine xstream) are checked later, once the engine is known
check_requirements(['pandas', 'Bio', 'metapredict', 'cython', 'matplotlib', 'protfasta'], not args.interactive)

## import modules from python3 standard library
import subprocess
import glob
from os import listdir
from os.path import isfile, join

## change import directory to this script's scripts/ folder for import, then import shared FLIPPer modules
sys.path.insert(1, os.path.join(SCRIPT_DIR, "scripts"))
from FLIPPer_lib import *

## captured before any chdir below, so --characterize (a standalone action resolved against the
## directory the user invoked FLIPPer from) is unaffected by --input-dir's chdir
INVOCATION_DIR = os.getcwd()

## move into the requested input directory (defaults to the current directory) so that all
## file discovery and output below operates on the data the user actually wants processed
os.chdir(args.input_dir)

## define path as current working directory
PATH = os.getcwd()

## Find all non-directory files in the current directory, then remove FLIPPer default files and OS-specfic hidden files from list
onlyfiles = [f for f in listdir(PATH) if isfile (join(PATH,f))]
package_files = ["FLIPPer.py", "README.md", "Changelog.txt", "requirements.txt", "desktop.ini", ".DS_Store", ".gitattributes", "LICENSE", "icon.svg", ".gitignore"]
for file in package_files:
    if os.path.isfile(file):
            onlyfiles.remove(file)

## Check if output folder for files exists
## If it does, exclude from files for analysis then continue
## (built as a new list rather than mutating onlyfiles while iterating over it,
## which used to silently skip files when two in a row already had output folders)
already_processed = [file for file in onlyfiles if os.path.exists("{}_FLIPPer_outputs".format(file))]
for file in already_processed:
    print(lineenter)
    print("Output folder for "+str(file)+" already exists, please rename or remove - skipping.")
onlyfiles = [f for f in onlyfiles if f not in already_processed]

## Print blurb
print(lineenter)
print("FLIPPer - Fast Linker Identification Pipeline for Pyrenoids - v3.0 (J. Barrett [james.barrett@york.ac.uk])")
print(lineenter)

## Which repeat-detection engine to use - from --engine if given, defaulted to xstream if
## omitted (preserves old scripts/HPC jobs that predate --engine), otherwise, under --interactive,
## asked here before any other prompt or dependency check below depends on knowing it
if args.engine:
    engine_name = args.engine
elif args.interactive:
    engine_choice = input("Which repeat-detection engine? (xstream/detectrepeats) [xstream]: ").strip().lower()
    engine_name = engine_choice if engine_choice in ("xstream", "detectrepeats") else "xstream"
else:
    engine_name = "xstream"

if engine_name == "xstream":
    import engine_xstream as engine
    from engine_xstream import *
else:
    import engine_detectrepeats as engine
    from engine_detectrepeats import *

## this engine's own extra Python package requirements (e.g. bs4 for xstream)
check_requirements(engine.EXTRA_PY_REQUIREMENTS, not args.interactive)

## Check that the engine's own external dependency (Java for XSTREAM, Rscript+DECIPHER for
## DetectRepeats) is actually available before asking the user anything else - every file's
## processing depends on it, and failing fast avoids wasting the user's time on prompts before a
## guaranteed failure.
if not engine.check_dependencies():
    print("Exiting.")
    sys.exit()

## --characterize is a standalone action: characterize the given target/reference FASTA and exit,
## without touching --input-dir or running the main proteome-scanning pipeline below
if args.characterize:
    target_path = args.characterize if os.path.isabs(args.characterize) else os.path.join(INVOCATION_DIR, args.characterize)
    target_path = os.path.abspath(target_path)
    if not os.path.isfile(target_path):
        print(target_path + " not found. Exiting.")
        sys.exit(1)
    if not validate_fasta(target_path):
        print(target_path + " is not FASTA format. Exiting.")
        sys.exit(1)
    os.chdir(os.path.dirname(target_path))
    engine.characterize(os.path.basename(target_path))
    sys.exit(0)

## --refilter is a standalone action, like --characterize: re-apply post-detection filters to an
## existing DetectRepeats output folder and exit, without touching --input-dir or running the main
## proteome-scanning pipeline below. Only detectrepeats saves the raw, pre-filter report this reads
## (see engine_detectrepeats.py's process_file/refilter) - xstream has no equivalent yet.
if args.refilter:
    if engine_name != "detectrepeats":
        print("--refilter is only supported with --engine detectrepeats (pass --engine detectrepeats). Exiting.")
        sys.exit(1)
    refilter_dir = args.refilter if os.path.isabs(args.refilter) else os.path.join(INVOCATION_DIR, args.refilter)
    refilter_dir = os.path.abspath(refilter_dir)
    if not os.path.isdir(refilter_dir):
        print(refilter_dir + " not found. Exiting.")
        sys.exit(1)
    metapredict_plot_flag = 'y' if args.metapredict_plot else 'n'
    ok = engine.refilter(refilter_dir, metapredict_plot_flag, args.metapredict_filter_value, args.aromatic, args.electrostatic,
                          args.min_copies, args.min_period, args.max_period, args.coverage)
    sys.exit(0 if ok else 1)

if args.interactive:
    ## Interactive wizard: one line per parameter, showing the flag/default value as the bracketed
    ## default - press Enter to keep it, or type a replacement. This replaces the old two-step
    ## "customize this group? y/n" gate followed by a separate list of inputs.
    def ask(prompt_text, default, cast=str):
        raw = input("{} [{}]: ".format(prompt_text, default)).strip()
        return cast(raw) if raw else default

    print(lineenter + "\n" + "Set parameters - press Enter to keep the default shown in brackets..." + "\n")

    pI = ask("pI threshold", args.pi, float)
    pI_direction = ask("pI direction - keep pI >= threshold ('min', for basic linkers like EPYC1/CsLinker) or pI <= threshold ('max', for acidic targets) (min/max)", args.pi_direction)
    if pI_direction not in ("min", "max"):
        pI_direction = args.pi_direction
    THRatio = ask("Turn/Helix ratio threshold", args.th_ratio, float)
    Serine = ask("Serine content threshold", args.serine, float)
    Alanine = ask("Alanine content threshold", args.alanine, float)

    print(lineenter)
    if engine_name == "xstream":
        Copy = ask("XSTREAM minimum copy number", args.copy)
        Word = ask("XSTREAM minimum word match", args.word)
        Consensus = ask("XSTREAM consensus match", args.consensus)
        Gaps = ask("XSTREAM maximum gaps in repeats", args.gaps)
    else:
        MinScore = ask("DetectRepeats minimum repeat significance score", args.min_score)
        MinCopies = ask("Minimum tandem repeat copy number", args.min_copies)
    minPeriod = ask("Minimum repeat period (aa)", args.min_period)
    maxPeriod = ask("Maximum repeat period (aa)", args.max_period)
    Coverage = ask("Sequence proportion repeats must cover (0.4 for fusions, 0.75 for full linkers)", args.coverage)

    print(lineenter)
    Aromatic = ask("Minimum aromatic residues in repeat region", args.aromatic, float)
    Electrostatic = ask("Minimum electrostatic residues in repeat region", args.electrostatic, float)

    print(lineenter)
    full_output = ask("Keep full sequence analysis of all input sequences? - large file, n recommended (y/n)", "y" if args.full_output else "n")
    metapredict_filter_value = ask("Minimum metapredict disorder % required within the identified repeat region", args.metapredict_filter_value, float)
    metapredict_plot = ask("Plot metapredict/pLDDT profiles for candidate sequences? - y recommended (y/n)", "y" if args.metapredict_plot else "n")
else:
    ## Default: use the CLI-provided (or documented default) parameters directly, no prompts
    pI = args.pi
    pI_direction = args.pi_direction
    THRatio = args.th_ratio
    Serine = args.serine
    Alanine = args.alanine
    if engine_name == "xstream":
        Copy = args.copy
        Word = args.word
        Consensus = args.consensus
        Gaps = args.gaps
    else:
        MinScore = args.min_score
        MinCopies = args.min_copies
    minPeriod = args.min_period
    maxPeriod = args.max_period
    Coverage = args.coverage
    Aromatic = args.aromatic
    Electrostatic = args.electrostatic
    full_output = "y" if args.full_output else "n"
    metapredict_filter_value = args.metapredict_filter_value
    metapredict_plot = "y" if args.metapredict_plot else "n"

print("Running with the following parameters:")
print("\tEngine: ", engine_name)
print("\tpI threshold: ", pI, "(direction: {})".format(pI_direction))
print("\tTurn/Helix ratio threshold: ", THRatio)
print("\tSerine content threshold: ", Serine)
print("\tAlanine content threshold: ", Alanine)
if engine_name == "xstream":
    print("\tXSTREAM minimum copy number: ", Copy)
    print("\tXSTREAM minimum word match: ", Word)
    print("\tXSTREAM consensus match: ", Consensus)
    print("\tXSTREAM maximum gaps: ", Gaps)
else:
    print("\tDetectRepeats minimum score: ", MinScore)
    print("\tMinimum copy number: ", MinCopies)
print("\tMinimum repeat period: ", minPeriod)
print("\tMaximum repeat period: ", maxPeriod)
print("\tCoverage: ", Coverage)
print("\tMinimum aromatic residues: ", Aromatic)
print("\tMinimum electrostatic residues: ", Electrostatic)
print("\tKeep full sequence analysis output: ", full_output)
print("\tmetapredict filter value: ", metapredict_filter_value)
print("\tPlot metapredict/pLDDT profiles: ", metapredict_plot)
print(lineenter)

## for each file in onlyfiles without output directory already existing
for file in onlyfiles:
    print(file)
    ## validate if the file is fasta, if not, skip it
    if not validate_fasta(file):
        print(lineenter)
        print(file + " is not FASTA format - skipping.")
        print(lineenter)
        continue

    ## make directories for outputs
    destination_folder="{}_FLIPPer_outputs".format(file)
    directory = "{}_FLIPPer_outputs/metapredict_plots".format(file)
    os.mkdir(destination_folder)

    if metapredict_plot == 'y':
        os.makedirs(directory)

    ## Everything below is wrapped in try/finally: any failure (the engine erroring out, no
    ## repeats found, an unexpected exception) is reported and this file is skipped, but temp
    ## files are always cleaned up and whatever was produced is always moved into
    ## destination_folder - so a failure never leaves debris that corrupts a later run.
    try:
        ## run analysis and filtering module (shared, engine-independent)
        analysis_and_filtering(file, pI, THRatio, Serine, Alanine, Aromatic, Electrostatic, full_output, pI_direction)

        ## run the chosen engine's own search/extract/metapredict/report pipeline
        if engine_name == "xstream":
            reached_end = engine.process_file(file, PATH, directory, metapredict_plot, metapredict_filter_value, Aromatic, Electrostatic,
                                               Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage)
        else:
            reached_end = engine.process_file(file, PATH, directory, metapredict_plot, metapredict_filter_value, Aromatic, Electrostatic,
                                               MinScore, MinCopies, minPeriod, maxPeriod, Coverage)
        if not reached_end:
            continue

        print("Done!")
        print(lineenter)

        ## pass variables from input to output_variables module to write variables file
        if engine_name == "xstream":
            engine.output_variables(file, pI, THRatio, Serine, Alanine, Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage, Aromatic, Electrostatic, metapredict_filter_value, pI_direction)
        else:
            engine.output_variables(file, pI, THRatio, Serine, Alanine, MinScore, MinCopies, minPeriod, maxPeriod, Coverage, Aromatic, Electrostatic, metapredict_filter_value, pI_direction)

        ## end of analysis run
        print("Cleaning up files from "+str(file)+ " analysis.")
        print("Analysis of "+str(file)+" finished!")
        print(lineenter)
    except Exception as exc:
        print(lineenter)
        print("Unexpected error while processing " + str(file) + ": " + str(exc))
        print("Skipping to next file.")
        print(lineenter)
    finally:
        ## always move whatever was produced into destination_folder, then remove genuinely
        ## disposable scratch files - in that order, since cleanup_temp_files() also sweeps up
        ## stray engine temp files and must not run before finalize_output() has had a chance to
        ## save the ones that are real output (this run's final engine report)
        finalize_output(file, destination_folder, PATH, engine)
        cleanup_temp_files(engine)

print("ALL FINISHED!")
