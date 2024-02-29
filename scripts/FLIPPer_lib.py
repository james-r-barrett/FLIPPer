## Define graphics for print outputs
lineenter = "\n" + "=========================================================" + "\n" 

## Default variables
## Specify default values, these are overwritten by user input if y answered
pI= "8"
THRatio= "1"
Serine= "0.05"
Alanine= "0.01"
Copy= "3"
Word= "0.3625"
Consensus= "0.4"
Gaps= "50"
minPeriod= "20"
maxPeriod= "120"
Coverage= "0.75"
Type= "long"
Aromatic= "1"
Electrostatic= "2"

## module to check if each of the files in the directory are fasta formatted
## returns true if is fasta file
def validate_fasta(filename):
    from Bio import SeqIO
    with open(filename, "r", encoding='ascii') as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)

## module to complete analysis of input sequences with ProtParam from biopython
## then filter to input variables using pandas dataframe
def analysis_and_filtering(file, pI, THratio, Serine, Alanine, Aromatic, Electrostatic, full_output):
    import os
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import pandas as pd
    import re
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
    print("Once complete, the outputs will be placed in {}_outputs".format(file))
    for record in SeqIO.parse(file,'fasta'): #input the file from master script as "f", in fasta format
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
        AA_out.append(Y.get_amino_acids_percent()) #gives you the aminoacid percentage of all the aminoacids
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
    if not os.path.exists('Sequence_Analysis'):
        os.mkdir('Sequence_Analysis')
        print("Temporary directory " "'""Sequence_Analysis""'" " created")
    else:
        print("Directory ""'""Sequence_Analysis""'" "already exists, please analyse output carefully")
    ## Create filename from input by removing extension with regex command
    no_extension = re.match(r"(.*)\.",file).group(0)
    ## If user specifies to keep full sequence analysis, write to file
    if full_output == "y":
        df.to_csv('Sequence_Analysis/%s_FullAnalysis.txt' % no_extension, index = None, sep='\t')
    ## Create filtered variable table and write to file
    Ratiodf = df[df['Ratio (T/H)']>=float(THratio)]
    RatiodfSerine= Ratiodf.loc[(df['%S'])>=float(Serine)]
    RatiodfSerineAlanine= RatiodfSerine.loc[df['%A']>=float(Alanine)]
    RatiodfPI = RatiodfSerineAlanine.loc[df['pI']>=float(pI)]    
    RatiodfPI.to_csv('Sequence_Analysis/%s_FilteredCandidates.TSV' % no_extension, index=None, sep='	')
    ## Filter sequences and IDs from dataframe
    ## Then write to temporary file used for Xstream
    FASTA=RatiodfPI.iloc[:,0:2]
    FASTA.to_csv('CandidateSequences_Temp.FASTA', index=None, header=None, sep='\n')
    ## Get number of sequence in input and output
    ## Print to user numbers
    seqno = (len(df.index))
    filteredseqno = (len(RatiodfPI.index))
    print("Sequence analysis of " +str(seqno)+ " sequences complete!")
    print(str(filteredseqno) + " sequences after filtering, passing to Xstream.")
    print(lineenter)

## module to extract candidate sequences from html3 file output from Xstream
## then optionally filter to only include sequences with at least one aromatic and three electrostatic residues in repeat region
def xstream_extract2(f, input_file, Aromatic, Electrostatic):
    from bs4 import BeautifulSoup
    import sys
    import fileinput
    import re
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    ## lists for repeat regions and their IDs
    repeats=[]
    IDs=[]
    with open(f) as html_file:
        ## use bautifulsoup html parser to extract text from _2.html file passed from FLIPPer.py
        soup= BeautifulSoup(html_file, "html.parser")
        html_file.close()
        content=soup.get_text()
        ## use regex to extract the consensus repeat motif
        ## separately, use regex to extract ID
        repeat_sequence = re.findall("=\n[\D]*\n", content)
        ID = re.findall(".*Position", content)
        ## clean up regex extractions
        for item in repeat_sequence:
            rep_seq = re.sub("[:|\*| ]*\n", "", item)
            rep_seq2 = re.sub("=","", rep_seq)
            repeats.append(rep_seq2)
        ## clean up regex extraction of ID
        for item in ID:
            ID2 = re.sub("Position","", item)
            IDs.append(ID2)
        ## add identifier, then create pandas dataframe with ID and consensus repeat motif, then write to temp file in fasta format
        IDs_fasta = [">" + ID for ID in IDs]
        results = zip(IDs_fasta, repeats)
        df_xstream = pd.DataFrame(results)
        df_xstream.to_csv('Temp_Xstream_positives.fasta', index=None, header=None, sep='\n')
    ## lists for post-xstream filtering
    xstream_sequence =[]
    xstream_AA = []
    xstream_AA_val =[]
    xstream_ID=[]
    number_before_filterting=[]
    for record in SeqIO.parse("Temp_Xstream_positives.fasta",'fasta'): #input the file from master script as "f", in fasta format
        number_before_filterting.append(record)
        ## Import fasta format temp file with ID and consensus repeat outputted from above
        ## Then output sequence to pre-defined list
        xsequence=format(record.seq)
        xstream_sequence.append(xsequence)
        Y = ProteinAnalysis(xsequence)
        ## output counted amino acids for each motif
        xstream_AA.append(Y.count_amino_acids())
        xstream_ID.append(format(record.id))   
    ## extract just values from protparam output
    for listitem in xstream_AA:
        xstream_AA_val.append(list(listitem.values()))
    ## create dataframe and input AA count in labelled columns
    ## then insert ID extracted from xstream output
    ## then create column for aromatic and electrostatic count per consensus
    xstream_df = pd.DataFrame(xstream_AA_val, columns = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])
    xstream_df.insert(0, 'ID', xstream_ID)
    xstream_df.insert(1, 'Aromatics', (xstream_df['F']+xstream_df['W']+xstream_df['Y']))
    xstream_df.insert(2, 'Electrostatics', (xstream_df['D']+xstream_df['E']+xstream_df['R']+xstream_df['K']))
    ## filter datafram to threshold value (either default or inputted)
    xstream_dfaromatic = xstream_df[xstream_df['Aromatics']>=float(Aromatic)]
    xstream_dfaromaticelectrostatic = xstream_dfaromatic.loc[xstream_df['Electrostatics']>=float(Electrostatic)]
    ## output the IDs of the filtered proteins to a list
    positive_IDs = list(xstream_dfaromaticelectrostatic['ID'])
    ## creat new temp file that contains fasta format of IDs of filtered proteins with sequence extracted from input_file passed from FLIPPer.py
    number_after_filtering=[]
    with open ("Temp_xstream_filtered.fasta", "w") as f:
        for record in SeqIO.parse(input_file, "fasta"):
            if record.id in positive_IDs:
                SeqIO.write([record],f,"fasta")
                number_after_filtering.append(record)
    print(lineenter)
    print("Number of sequences before filtering repeat regions:",len(number_before_filterting))
    print("Number of sequences after filtering repeat regions:",len(number_after_filtering))
    print(lineenter)

def metapredict_htp(file_name, directory, metapredict_plot, metapredict_filter_value):
    import protfasta
    import re
    import metapredict as meta
    import pandas as pd
    protfasta_seqs = protfasta.read_fasta(file_name, invalid_sequence_action = 'convert', return_list = True)
    IDs = []
    sequences = []
    disorder_score = []
    for seqs in protfasta_seqs:
        IDs.append(seqs[0])
        sequences.append(seqs[1])
        disorder_score.append(meta.percent_disorder(seqs[1]))
    fasta_IDs = [">" + ID for ID in IDs]
    dict = {'IDs': fasta_IDs, 'seq': sequences, 'score': disorder_score} 
    df=pd.DataFrame(dict)
    if metapredict_plot == 'y':
        filtered_df=df[df['score']>50]
        filtered_cands=filtered_df.iloc[:,0:2]
        filtered_cands.to_csv('candidate_sequences.fasta', index=None, header=None, sep='\n')
        filtered_cands.to_csv('candidate_sequences.csv', index=None, header=None, sep=',')
        protfasta_seqs = protfasta.read_fasta("candidate_sequences.fasta", invalid_sequence_action = 'convert', return_list = True)
        n=len(protfasta_seqs)
        i=0
        print("Initialising metapredict module...")
        print('Number of sequences for metapredict filtering/plotting: '+str(n)+ '.\n')
        for seqs in protfasta_seqs:
            PlotID=re.sub("\||\*|\?|\.|\/|\"|\<|\>|\:","_",seqs[0])
            PlotID=PlotID.split(' ')[0]
            meta.graph_disorder(seqs[1], pLDDT_scores=True, DPI=300, output_file=directory+'/%s_metapredict_plot.png' %PlotID, title = "%s" %seqs[0])
            i += 1      
            print("Plotted " +str(i)+" sequences.")

##module to output variables
def output_variables(file, pI, THRatio, Serine, Alanine, Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage, Type, Aromatic, Electrostatic, metapredict_filter_value):
    import sys    
    print("Outputting variables used.")
    tem = sys.stdout
    sys.stdout= m =open('{}_variables.txt'.format(file),'w')
    print("==========================================")
    print("Filtering variables:")
    print("\tpI Threshold: ", pI)
    print("\tTurn/Helix Ratio Threshold: ", THRatio)
    print("\tSerine content threshold: ", Serine)
    print("\tAlanine content threshold: ", Alanine)
    print("==========================================")
    print("XSTREAM variables:")
    print("\tMinimum copy number: ", Copy)
    print("\tMinmum word match: ", Word)
    print("\tConsensus match: ", Consensus)
    print("\tMaximum gaps in repeats: ", Gaps)
    print("\tMinimum repeat period: ", minPeriod)
    print("\tMaximum repeat period: ", maxPeriod)
    print("\tSequence coverage: ", Coverage)
    print("==========================================")
    print("Post-XSTREAM Filtering")
    print("\tNumber of aromatic residues: ", Aromatic)
    print("\tNumber of electrostatic residues: ", Electrostatic)
    print("==========================================")
    print("metapredict Filtering")
    print("\tmetapredict filtering value: ", metapredict_filter_value)
    sys.stdout = tem
    m.close()

## module to move files after pipeline run       
def move_files(source_file_name, destination_folder_name, path):
    import os
    os.rename("{}/{}".format(path, source_file_name), "{}/{}/{}".format(path, destination_folder_name, source_file_name))
