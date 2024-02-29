=============================================================
Fast Linker Idenfitification Pipeline for Pyrenoids (FLIPPer)
=============================================================

FLIPPer runs against FASTA formatted sequences. The FASTA formatted files must be in the FLIPPer directory root directory (i.e. where the FLIPPer.py script is).
The pipeline will scan the root directory and check for files formatted in the FASTA format then run the pipeline against them using the defined settings.


Usage:

python3 FLIPPer.py

Runnning the above command from terminal or command prompt will lanuch the paramterisation interface.


Outputs:

For each of the analysed FASTA files several output are produced.

1. IUPred plots
	- Overlaid PNGs of disorder prediction, net charge per residue (NCPR) and hydrophobicity profiles for each entry in the FASTA file.
2. XSTREAM outputs
	- The identified tandem repeats in each of the sequences, viewable as a linked HTML file (suffix _out_2.html)
3. Candidate sequences
	- A FASTA file of all of the outputted candidate linker sequences.
	- A TSV file of the same sequences
4. Variables file
	- A text file containing the user-defined variables used in the pipeline run.
5. Sequences analysis
	- TSV files of the sequence analysis completed as step 1 of the sequence filtering


Requirements:

python 3.10 or higher
biopython
matplotlib
beautifulsoup4
localcider
pandas


This pipeline utilises XSTREAM, a tandem repeat detection software developed by Aaron Newman and James Cooper. If you make use of this pipeline, please be sure to also credit their work:

Newman, A. M. and Cooper, J. B. (2007). XSTREAM: a practical algorithm for identification and architecture modeling of tandem repeats in protein sequences. BMC bioinformatics, 8, p.382.
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-382

This pipeline utilises metapredict, an intrinsic disorder prediction software developed by Ryan Emenecker, Daniel Griffith and Alex Holehouse, please be sure to also credit their work:

Emenecker, R. J., Griffith, D. and Holehouse, A. S. (2022). Metapredict V2: An update to metapredict, a fast, accurate, and easy-to-use predictor of consensus disorder and structure. bioRxiv, p.2022.06.06.494887.
https://www.biorxiv.org/content/10.1101/2022.06.06.494887v1

Emenecker, R. J., Griffith, D. and Holehouse, A. S. (2021). Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure. Biophysical journal, 120 (20), pp.4312–4319.
https://www.cell.com/biophysj/fulltext/S0006-3495(21)00725-6

Note:

The record id (the FASTA name) must not exceed 50 characters, as this will interfere with sequence filtering after XSTREAM repeat detection.
