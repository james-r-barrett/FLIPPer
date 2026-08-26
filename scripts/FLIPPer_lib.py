## re-exported so `from FLIPPer_lib import *` (as FLIPPer.py does) also picks up the report
## builder shared with both engines - see scripts/candidate_report.py
from candidate_report import build_candidate_report

## Define graphics for print outputs
lineenter = "\n" + "=========================================================" + "\n"

## Default variables shared by both repeat-detection engines
## Specify default values, these are overwritten by user input if y answered
pI= "8"
pI_direction = "min"
THRatio= "1"
Serine= "0.05"
Alanine= "0.01"
minPeriod= "20"
maxPeriod= "120"
Coverage= "0.75"
Aromatic= "1"
Electrostatic= "2"

## module to check if each of the files in the directory are fasta formatted
## returns true if is fasta file
def validate_fasta(filename):
    from Bio import SeqIO
    try:
        with open(filename, "r", encoding='utf-8') as handle:
            fasta = SeqIO.parse(handle, "fasta-pearson")
            return any(fasta)
    except (UnicodeDecodeError, ValueError):
        return False

## module to complete analysis of input sequences with ProtParam from biopython
## then filter to input variables using pandas dataframe
def analysis_and_filtering(file, pI, THratio, Serine, Alanine, Aromatic, Electrostatic, full_output, pI_direction="min"):
    import os
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import pandas as pd
    import fileinput
    print(lineenter)
    ## Define empty lists to use in below analysis
    seq_out=[]
    length_out=[]
    AA_out=[]
    pI_out=[]
    ID_out=[]
    values_out=[]
    print("Running sequence analysis of " "'"+file+"'"+"\n")
    print("Temporary files will be created within the directory, please do not remove them.")
    print("Once complete, the outputs will be placed in {}_FLIPPer_outputs".format(file))
    for record in SeqIO.parse(file,'fasta-pearson'): #input the file from master script as "f", in fasta format
        ## Import sequences from "file" passed from master script as f
        ## Then output sequecne to pre-defined list
        ## And output length of sequence to separate list
        sequence=format(record.seq)
        seq_out.append(sequence)
        length_out.append(len(sequence))
        ## Use ProteinAnalysis module from Bio.SeqUtils.ProtParam (biopython dependency) to get sequnce properties
        ## Define sequence from above as string
        ## Then get percentage of each amino acid
        ## Then get the isoelectric point of the sequence
        ## And output ID to separate list
        Y = ProteinAnalysis(sequence)
        ## Y.amino_acids_percent is on a 0-100 scale in current Biopython versions, but the
        ## Serine/Alanine thresholds below (and the Aromatic/Electrostatic residue-count columns
        ## further down, which multiply these by Length) assume a 0-1 fraction - convert here so
        ## --serine/--alanine keep meaning "fraction of residues", not "percent >= 0.05"
        AA_out.append({aa: value / 100.0 for aa, value in Y.amino_acids_percent.items()}) #gives you the aminoacid fraction of all the aminoacids
        pI_out.append(Y.isoelectric_point())
        ID_out.append('>' + format(record.id))
    ## Extract values from AA_out list, which contain data in format A: XXXX, where XXXX is desired info
    for listitem in AA_out:
        values_out.append(list(listitem.values()))
    ## Create data frame from values, with column titles for each amino acid
    df = pd.DataFrame(values_out, columns = ['%A','%C','%D','%E','%F','%G','%H','%I','%K','%L','%M','%N','%P','%Q','%R','%S','%T','%V','%W','%Y'])
    ## Insert the rest of the desired data, also creating columns with added/divded values where required - from ExPasy?
    df.insert(0,'ID', ID_out)
    df.insert(1,'Sequence', seq_out)
    df.insert(2,'Length', length_out)
    df.insert(3,'Aromatic',((df['%F']*df['Length'])+(df['%W']*df['Length'])+(df['%Y']*df['Length'])))
    df.insert(4,'Electrostatic',((df['%R']*df['Length'])+(df['%K']*df['Length'])+(df['%D']*df['Length'])+(df['%E']*df['Length'])))
    df.insert(5,'Fraction Expanding', (df['%R'] + df['%K'] + df['%D'] + df['%E'] + df['%P']))
    df.insert(6,'Fraction Disorder Promoting', (df['%A'] + df['%G'] + df['%R'] + df['%D'] + df['%H'] + df['%Q'] + df['%K'] + df['%S'] + df['%E'] + df['%P']))
    df.insert(7,'Helix%', (df['%F'] + df['%I'] + df['%L'] + df['%V'] + df['%W'] + df['%Y']))
    df.insert(8,'Turn%', df['%P'] + df['%N'] + df['%G'] + df['%S'])
    df.insert(9,'Sheet%', df['%E'] + df['%M'] + df['%A'] + df['%L'])
    df.insert(10,'Ratio (T/H)',(df['Turn%']/df['Helix%']))
    df.insert(11,'Charged Ratio (Positive/Negative)', ((df['%R'] + df['%K'])/(df['%D'] + df['%E'])))
    df.insert(12,'pI', pI_out)
    ## Create Sequence Analysis directory to put output in
    if not os.path.exists('sequence_analysis'):
        os.mkdir('sequence_analysis')
        print("Temporary directory " "'""sequence_analysis""'" " created")
    else:
        print("Directory ""'""sequence_analysis""'" "already exists, please analyse output carefully")
    ## Create filename from input by removing extension
    no_extension = os.path.splitext(file)[0]
    ## If user specifies to keep full sequence analysis, write to file
    if full_output == "y":
        df.to_csv('sequence_analysis/%s_FullAnalysis.txt' % no_extension, index = None, sep='\t')
    ## Create filtered variable table and write to file
    Ratiodf = df[df['Ratio (T/H)']>=float(THratio)]
    RatiodfSerine= Ratiodf.loc[(df['%S'])>=float(Serine)]
    RatiodfSerineAlanine= RatiodfSerine.loc[df['%A']>=float(Alanine)]
    ## pI_direction="min" (default) keeps basic proteins (pI >= threshold) - the profile of
    ## known pyrenoid linkers like EPYC1/CsLinker. pI_direction="max" flips this to keep acidic
    ## proteins (pI <= threshold) instead, for target families with the opposite composition.
    if str(pI_direction).lower() == "max":
        RatiodfPI = RatiodfSerineAlanine.loc[df['pI']<=float(pI)]
    else:
        RatiodfPI = RatiodfSerineAlanine.loc[df['pI']>=float(pI)]
    RatiodfPI.to_csv('sequence_analysis/%s_FilteredCandidates.csv' % no_extension, index=None, sep=',')
    ## Filter sequences and IDs from dataframe
    ## Then write to temporary file used for the repeat-detection engine
    FASTA = RatiodfPI.iloc[:, 0:2]
    # Write to FASTA file
    with open('CandidateSequences_Temp.FASTA', 'w') as f:
        for i, row in FASTA.iterrows():
            f.write(f">{row['ID']}\n{row['Sequence']}\n")
    ## Get number of sequence in input and output
    ## Print to user numbers
    seqno = (len(df.index))
    filteredseqno = (len(RatiodfPI.index))
    print("Sequence analysis of " +str(seqno)+ " sequences complete!")
    print(str(filteredseqno) + " sequences after filtering, passing to repeat detection.")
    print(lineenter)

## module to write out the final candidate set (already filtered upstream, per-repeat-region, on
## composition and metapredict disorder - see xstream_extract2/detect_repeats_extract) to
## candidate_sequences.fasta/.csv, and optionally plot each one's whole-protein disorder/pLDDT
## profile for visual inspection - engine-independent, shared by both engines
def metapredict_htp(file_name, directory, metapredict_plot):
    import protfasta
    import re
    import metapredict as meta
    import pandas as pd
    protfasta_seqs = protfasta.read_fasta(file_name, invalid_sequence_action = 'convert', return_list = True)
    IDs = []
    sequences = []
    for seqs in protfasta_seqs:
        IDs.append(seqs[0])
        sequences.append(seqs[1])
    fasta_IDs = [">" + ID for ID in IDs]
    dict = {'IDs': fasta_IDs, 'seq': sequences}
    df=pd.DataFrame(dict)
    with open('candidate_sequences.fasta', 'w') as f:
        for _, row in df.iterrows():
            f.write(f"{row['IDs'].strip()}\n{row['seq'].strip()}\n")
    df.to_csv('candidate_sequences.csv', index=None, header=None, sep=',')
    candidate_number = len(df.index)
    print(str(candidate_number)+' final candidates written to candidate_sequences.fasta/csv.\n')
    if metapredict_plot == 'y':
        protfasta_seqs = protfasta.read_fasta("candidate_sequences.fasta", invalid_sequence_action = 'convert', return_list = True)
        i=0
        print('Number of sequences for metapredict filtering/plotting: '+str(candidate_number)+ '.\n')
        for seqs in protfasta_seqs:
            PlotID = re.sub(r"[|*?./\"<>:]", "_", seqs[0])
            PlotID=PlotID.split(' ')[0]
            meta.graph_disorder(seqs[1], pLDDT_scores=True, DPI=300, output_file=directory+'/%s_metapredict_plot.pdf' %PlotID, title = "%s" %seqs[0])
            i += 1
            print("Plotted " +str(i)+" sequences.")

## module to move files after pipeline run
def move_files(source_file_name, destination_folder_name, path):
    import os
    os.rename("{}/{}".format(path, source_file_name), "{}/{}/{}".format(path, destination_folder_name, source_file_name))

## module to remove known temporary files left over from a (possibly failed) pipeline run
## so that they can never be picked up as input files on a later run of FLIPPer.py - the shared
## temp file plus whatever engine-specific ones the chosen engine module knows about
def cleanup_temp_files(engine):
    import os
    if os.path.exists("CandidateSequences_Temp.FASTA"):
        os.remove("CandidateSequences_Temp.FASTA")
    engine.cleanup_extra()

## module to move whatever output artifacts exist for a file into its destination folder
## called unconditionally at the end of each file's processing (success or failure) so that
## partial results never linger in, or pollute, the working directory - common artifacts (both
## engines produce these) plus whatever engine-specific ones the chosen engine module knows about
def finalize_output(file, destination_folder, path, engine):
    import os
    import glob
    engine.finalize_extra(destination_folder, path)
    for z in glob.glob("*_candidate_report.html"):
        move_files(z, destination_folder, path)
    if os.path.exists("sequence_analysis"):
        move_files("sequence_analysis", destination_folder, path)
    if os.path.exists(file):
        move_files(file, destination_folder, path)
    if os.path.exists(file+"_variables.txt"):
        move_files(file+"_variables.txt", destination_folder, path)
    if os.path.exists("candidate_sequences.fasta"):
        move_files("candidate_sequences.fasta", destination_folder, path)
    if os.path.exists("candidate_sequences.csv"):
        move_files("candidate_sequences.csv", destination_folder, path)
