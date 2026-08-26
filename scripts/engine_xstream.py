## XSTREAM repeat-detection engine for FLIPPer.py - see engine_detectrepeats.py for the other
## engine. Both modules expose the same small interface (NAME, EXTRA_PY_REQUIREMENTS,
## check_dependencies, process_file, characterize, output_variables, finalize_extra,
## cleanup_extra) so FLIPPer.py's shared driver can call either interchangeably via --engine.
import os

NAME = "xstream"

## bs4 (BeautifulSoup) is only needed to parse XSTREAM's HTML reports - the DetectRepeats engine
## has no such extra requirement
EXTRA_PY_REQUIREMENTS = ['bs4']

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
XSTREAM_JAR = os.path.join(SCRIPT_DIR, "xstream.jar")

## Default variables for this engine's own search parameters, overwritten by user input if y
## answered - used both directly (non-interactive mode with no override) and, via
## `from engine_xstream import *`, as the interactive prompts' unmodified fallback
Copy= "3"
Word= "0.3625"
Consensus= "0.4"
Gaps= "55"

## checked once, right after the engine is chosen, before asking the user anything else - every
## file's processing depends on Java, and failing fast avoids wasting the user's time on prompts
## before a guaranteed failure
def check_dependencies():
    import subprocess
    try:
        java_ok = subprocess.call(["java", "-version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
    except FileNotFoundError:
        java_ok = False
    if not java_ok:
        print("Java could not be found on this system - XSTREAM requires a working Java installation.")
    return java_ok

## module to run XSTREAM as a subprocess, checking that it actually succeeded rather than assuming it did
## returns True on success, prints a message and returns False on failure (non-zero exit code or missing java)
def run_xstream(args, **kwargs):
    import subprocess
    try:
        result = subprocess.call(args, **kwargs)
    except FileNotFoundError as e:
        print("Could not run XSTREAM: " + str(e))
        return False
    if result != 0:
        print("XSTREAM exited with an error (code {}).".format(result))
        return False
    return True

## module to check whether an XSTREAM "_2.html" report actually contains any detected repeats.
## XSTREAM always writes this file, even when it finds nothing (the file just says literal text
## "No Repeats Found!") - so checking only that the file exists doesn't distinguish a genuine
## "nothing found" result from one with real candidates. Treating existence alone as "found"
## silently cascades empty temp files through the rest of the pipeline, and can crash the final
## XSTREAM report run later on (XSTREAM errors out on an empty FASTA input).
def xstream_found_repeats(html_files):
    for f in html_files:
        with open(f, encoding='utf-8', errors='ignore') as handle:
            content = handle.read()
        if "No Repeats Found" not in content:
            return True
    return False

## module to extract candidate sequences from html3 file output from Xstream, then filter to only
## include sequences with at least one aromatic and three electrostatic residues in the extracted
## consensus repeat motif AND at least metapredict_filter_value% disorder across the FULL repeat
## region the same repeat instance covers (its Begin-End span, all copies - not just the short
## consensus unit, which is a display convenience, not the region itself). Both are checked
## jointly against that same repeat instance, so a sequence with several repeat regions qualifies
## if one of them satisfies both together, not if different regions each satisfy a different
## criterion.
def xstream_extract2(f, input_file, Aromatic, Electrostatic, metapredict_filter_value):
    from bs4 import BeautifulSoup
    import sys
    import fileinput
    import re
    import os
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import metapredict as meta
    from protfasta.utilities import convert_to_valid
    from FLIPPer_lib import lineenter
    html_path = f  ## saved before `f` gets reused below as a temp-file handle variable
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
        repeat_sequence = re.findall(r"=\n[\D]*\n", content)
        ID = re.findall(".*Position", content)
        ## clean up regex extractions
        for item in repeat_sequence:
            rep_seq = re.sub(r"[:|\*| ]*\n", "", item)
            rep_seq2 = re.sub("=","", rep_seq)
            repeats.append(rep_seq2)
        ## clean up regex extraction of ID
        for item in ID:
            ID2 = re.sub("Position","", item)
            IDs.append(ID2)
        ## add identifier, then create pandas dataframe with ID and consensus repeat motif, then write to temp file in fasta format
        IDs_fasta = [">" + ID.lstrip(">") for ID in IDs]
        results = zip(IDs_fasta, repeats)
        with open('Temp_Xstream_positives.fasta', 'w') as f:
            for ID_line, repeat_seq in zip(IDs_fasta, repeats):
                f.write(f"{ID_line.strip()}\n{repeat_seq.strip()}\n")
    ## look up the full repeat-array span (Begin-End covering every copy, from XSTREAM's own
    ## Positions column) for each detected repeat, to score disorder over the actual repeat
    ## region rather than the short consensus unit above. single_sequence_fallback_id mirrors
    ## parse_xstream_repeats' handling of XSTREAM's single-sequence report layout, which omits
    ## the anchor parse_xstream_index otherwise reads the ID from - determined from whichever
    ## fasta this html report was actually generated against (the fixed temp filename XSTREAM was
    ## run against in process_file, still present on disk at this point in the pipeline).
    candidate_records = list(SeqIO.parse("CandidateSequences_Temp.FASTA", "fasta-pearson"))
    single_id = candidate_records[0].id if len(candidate_records) == 1 else None
    index_rows = parse_xstream_index(html_path, single_id)
    full_sequences = {r.id: str(r.seq) for r in SeqIO.parse(input_file, "fasta-pearson")}

    ## lists for post-xstream filtering
    xstream_sequence =[]
    xstream_AA = []
    xstream_AA_val =[]
    xstream_ID=[]
    xstream_disorder=[]
    number_before_filterting=[]
    for i, record in enumerate(SeqIO.parse("Temp_Xstream_positives.fasta",'fasta-pearson')): #input the file from master script as "f", in fasta format
        number_before_filterting.append(record)
        ## Import fasta format temp file with ID and consensus repeat outputted from above
        ## Then output sequence to pre-defined list
        xsequence=format(record.seq)
        xstream_sequence.append(xsequence)
        Y = ProteinAnalysis(xsequence)
        ## output counted amino acids for each motif
        xstream_AA.append(Y.count_amino_acids())
        xstream_ID.append(format(record.id))
        ## percent disorder of the full repeat region this same repeat instance covers (its
        ## Begin-End span) - index_rows comes from the same html report, in the same document
        ## order as this consensus-block extraction, so position i should always line up; the
        ## seq_id check guards against the two parses ever drifting out of step, degrading to the
        ## consensus unit itself for that one entry rather than misattributing another repeat's span
        span = index_rows[i] if i < len(index_rows) and index_rows[i]["seq_id"] == record.id else None
        if span is not None and record.id in full_sequences:
            region = full_sequences[record.id][int(span["begin"]) - 1:int(span["end"])]
        else:
            region = xsequence
        ## '-' gap characters (from XSTREAM's alignment display) and ambiguous residues (X/U/B/Z,
        ## common in real proteomes) are harmless for ProteinAnalysis.count_amino_acids above,
        ## which just ignores them, but metapredict validates its input strictly and raises on
        ## anything but the 20 standard amino acids - sanitize with the same conversion protfasta's
        ## invalid_sequence_action='convert' applies elsewhere in the pipeline (metapredict_htp)
        region_clean = convert_to_valid(region.upper())
        xstream_disorder.append(meta.percent_disorder(region_clean) if region_clean else 0.0)
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
    xstream_df.insert(3, 'Disorder', xstream_disorder)
    ## filter dataframe to threshold values (either default or inputted) - all three checked
    ## jointly against the same repeat instance/row
    xstream_dfaromatic = xstream_df[xstream_df['Aromatics']>=float(Aromatic)]
    xstream_dfaromaticelectrostatic = xstream_dfaromatic.loc[xstream_df['Electrostatics']>=float(Electrostatic)]
    xstream_dffiltered = xstream_dfaromaticelectrostatic.loc[xstream_df['Disorder']>float(metapredict_filter_value)]
    ## output the IDs of the filtered proteins to a list
    positive_IDs = list(xstream_dffiltered['ID'])
    ## creat new temp file that contains fasta format of IDs of filtered proteins with sequence extracted from input_file passed from FLIPPer.py
    number_after_filtering=[]
    with open ("Temp_xstream_filtered.fasta", "w") as f:
        for record in SeqIO.parse(input_file, "fasta-pearson"):
            if record.id in positive_IDs:
                SeqIO.write([record],f,"fasta")
                number_after_filtering.append(record)
    print(lineenter)
    print("Number of sequences before filtering repeat regions:",len(number_before_filterting))
    print("Number of sequences with a repeat region meeting "+str(Aromatic)+ ' aromatic (W/Y/F), ' +str(Electrostatic)+
          ' electrostatic (D/E/R/K) residues (in the consensus repeat unit), and >'+str(metapredict_filter_value)+
          '% disorder (across the full repeat region):',len(number_after_filtering))
    print(lineenter)

## module to extract the summary table (Positions/Period/CopyNumber/ConsensusError) from an
## XSTREAM "_2.html" report, in document order. XSTREAM lays this report out differently
## depending on the number of input sequences: with 2+ sequences each table is preceded by an
## <a href> link naming the sequence and has 4 data columns (Positions/Period/CopyNumber/
## ConsensusError); with exactly 1 sequence there's no such link (nothing to link to) and a 5th
## BlockColor column is added. Rather than branch on that, find the Positions cell by its
## "NN-NN" shape wherever it falls, and take the three cells that follow it as Period/CopyNumber/
## ConsensusError - this works for both layouts, including the header/data-row merging bug in
## the anchored one. single_sequence_fallback_id supplies the one sequence's ID for the
## no-anchor layout (pass None to skip rows FLIPPer can't attribute to a sequence).
def parse_xstream_index(html_path, single_sequence_fallback_id=None):
    from bs4 import BeautifulSoup
    import re
    position_pattern = re.compile(r'^\d+-\d+$')
    with open(html_path) as html_file:
        soup = BeautifulSoup(html_file, "html.parser")
    rows = []
    for table in soup.find_all('table', class_='MyFormat'):
        anchor = table.find_previous('a')
        ## a single-sequence report's nearest preceding <a> is a nameless jump target
        ## (<a name="1">), with no ID text of its own - anchor.get_text() is '' there, not None
        seq_id = anchor.get_text(strip=True) if anchor is not None else ""
        if seq_id:
            ## strip the "(length) ;type; coords" description XSTREAM appends after the ID,
            ## matching Biopython's record.id (just the first whitespace-delimited token)
            seq_id = seq_id.split()[0]
        else:
            if single_sequence_fallback_id is None:
                continue
            seq_id = single_sequence_fallback_id
        cells = [td.get_text(strip=True) for td in table.find_all('td')]
        pos_index = next((i for i, c in enumerate(cells) if position_pattern.match(c)), None)
        if pos_index is None or pos_index + 2 >= len(cells):
            continue
        positions, period, copy_number = cells[pos_index], cells[pos_index + 1], cells[pos_index + 2]
        consensus_error_cell = cells[pos_index + 3] if pos_index + 3 < len(cells) else None
        try:
            start, end = positions.split('-')
            start, end, period_val, copy_val = float(start), float(end), float(period), float(copy_number)
            error_val = float(consensus_error_cell) if consensus_error_cell not in (None, "") else None
        except ValueError:
            continue
        rows.append({
            "seq_id": seq_id, "begin": start, "end": end,
            "period": period_val, "copies": copy_val, "consensus_error": error_val,
        })
    return rows

## splits an XSTREAM "_2.html" report's raw HTML (not its parsed text - the alignment block
## below needs the literal per-residue <FONT> markup) into one chunk per "Repeat N" heading, in
## document order - the same order parse_xstream_index's rows come back in, so the two can be
## zipped together.
def _split_xstream_repeat_blocks(raw_html):
    import re
    starts = [m.start() for m in re.finditer(r'<b><FONT COLOR="FF0000">Repeat \d+', raw_html)]
    return [raw_html[s:e] for s, e in zip(starts, starts[1:] + [len(raw_html)])]

## pulls the per-residue character stream out of a run of XSTREAM's per-letter markup
## (<b><FONT COLOR="......">X</FONT></b>, one tag pair per residue/gap column, no whitespace
## separators) - this is how XSTREAM's report itself encodes both the full repeat-region
## alignment and the single-period reference row below it.
def _xstream_chars(segment):
    import re
    return re.findall(r'<b><FONT COLOR="[0-9A-Fa-f]{6}">(.)</FONT></b>', segment)

## Within one "Repeat N" block, XSTREAM shows the full repeat region as one long gapped
## alignment (its non-gap residue count matches End-Begin+1 for a multi-sequence report, though
## empirically it can run 1 residue over for the single-sequence report layout - unexplained, and
## clamped to `end` below rather than chased further), a "====" divider, then a single
## period-width reference/consensus row aligned to the same gap pattern. Empirically, len(reference
## row) * CopyNumber ~= len(full region row) - i.e. the full-region row is just the individual
## copies concatenated at reference-row width - so chunking it into reference-row-sized pieces
## reconstructs the same per-copy-aligned-against-consensus view DetectRepeats' alignment FASTAs
## show, including each copy's real sequence position (recovered by walking non-gap chars from
## Begin).
def _parse_xstream_alignment_block(segment, begin, end):
    import re
    divider_m = re.search(r'<FONT COLOR="0000FF">=+</FONT>', segment)
    if not divider_m:
        return [], [], [], ""
    pre_chars = _xstream_chars(segment[:divider_m.start()])
    rest = segment[divider_m.end():]
    match_m = re.search(r'<FONT COLOR="990000">', rest)
    ref_chars = _xstream_chars(rest[:match_m.start()] if match_m else rest)

    width = len(ref_chars)
    if width == 0 or not pre_chars:
        return [], [], [], ""

    copies, lefts, rights = [], [], []
    pos = int(begin)
    end = int(end)
    for i in range(0, len(pre_chars), width):
        chunk = pre_chars[i:i + width]
        non_gap = sum(1 for c in chunk if c != '-')
        if non_gap == 0:
            continue
        lefts.append(min(pos, end))
        pos = min(pos + non_gap, end + 1)
        rights.append(pos - 1)
        copies.append("".join(chunk) + "-" * (width - len(chunk)))
    return copies, lefts, rights, "".join(ref_chars)

## XSTREAM's ConsensusError (0-1, lower is better) is inverted into a 0-100 "match quality" score
## (higher is better) by parse_xstream_repeats below, so it sorts/colours through
## scripts/candidate_report.py's higher-is-better pill logic the same way DetectRepeats' score
## does - this is the tooltip text explaining that inversion on the report itself.
XSTREAM_SCORE_EXPLANATION = (
    "Match quality (100 x (1 - ConsensusError)): XSTREAM's own measure of how well this "
    "region's repeat copies agree with the single consensus repeat unit it fit to them, "
    "inverted here so higher = better, like the pill colouring elsewhere in this report. "
    "ConsensusError itself is 0 (perfect match) to 1."
)

## module to parse an XSTREAM "_2.html" report into the same normalized shape DetectRepeats'
## report CSV already has (see scripts/detect_repeats.R) - one row per detected repeat region,
## plus a sequence lookup and a matching alignment_provider - so both engines can render an
## identical candidate report through scripts/candidate_report.py's build_candidate_report() from
## nothing more than a couple of file paths.
def parse_xstream_repeats(html_path, candidate_fasta):
    import pandas as pd
    from Bio import SeqIO

    sequences = {r.id: str(r.seq) for r in SeqIO.parse(candidate_fasta, "fasta-pearson")}

    with open(html_path, encoding='utf-8', errors='ignore') as f:
        raw_html = f.read()

    single_id = next(iter(sequences)) if len(sequences) == 1 else None
    index_rows = parse_xstream_index(html_path, single_id)
    segments = _split_xstream_repeat_blocks(raw_html)
    ## these come from the same document in the same order, so should always match 1:1 - but if
    ## something about a report's layout isn't as expected, degrade to whichever prefix both
    ## agree on rather than risk misattributing one repeat's alignment to another's stats
    n = min(len(segments), len(index_rows))

    per_id_count = {}
    rows = []
    alignments = {}
    for info, segment in zip(index_rows[:n], segments[:n]):
        seq_id = info["seq_id"]
        seq = sequences.get(seq_id)
        if seq is None:
            continue
        repeat_index = per_id_count.get(seq_id, 0) + 1
        per_id_count[seq_id] = repeat_index

        copies, lefts, rights, consensus = _parse_xstream_alignment_block(segment, info["begin"], info["end"])
        aligned_records = [("copy{}".format(i + 1), s) for i, s in enumerate(copies)]
        if consensus:
            aligned_records.append(("consensus", consensus))
        alignments[(seq_id, repeat_index)] = aligned_records

        error_val = info["consensus_error"]
        rows.append({
            "ID": seq_id, "Begin": info["begin"], "End": info["end"],
            "Period": info["period"], "Copies": info["copies"],
            "Score": 100.0 * (1.0 - error_val) if error_val is not None else 0.0,
            "RepeatIndex": repeat_index,
            "Coverage": (info["end"] - info["begin"] + 1) / len(seq),
            "UnitLefts": ";".join(str(l) for l in lefts),
            "UnitRights": ";".join(str(r) for r in rights),
        })

    report_df = pd.DataFrame(rows, columns=[
        "ID", "Begin", "End", "Period", "Copies", "Score", "RepeatIndex", "Coverage", "UnitLefts", "UnitRights"])

    def alignment_provider(seq_id, repeat_index):
        return alignments.get((seq_id, repeat_index))

    return report_df, sequences, alignment_provider

## the per-file "core": run XSTREAM, extract/filter candidates, run metapredict, then rebuild the
## final XSTREAM report + candidate_report.html against the actual final candidate set. Called
## from FLIPPer.py's shared per-file loop, right after analysis_and_filtering(). Returns True if
## processing reached the end (whether or not a final report could be built), False if it bailed
## out early (XSTREAM run failed, or found no repeats at all) - the caller uses this the same way
## the old per-file loop used `continue`, to decide whether to print "Done!"/write the variables file.
def process_file(file, PATH, directory, metapredict_plot, metapredict_filter_value, Aromatic, Electrostatic,
                  Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage):
    import os
    import glob
    import subprocess
    from FLIPPer_lib import lineenter, validate_fasta, metapredict_htp, build_candidate_report

    ## run XSTREAM using filtered sequences
    if not run_xstream(["java", "-Xmx1000m", "-Xms1000m", "-jar", XSTREAM_JAR, "CandidateSequences_Temp.FASTA", "-e"+str(Copy), "-i"+str(Word), "-I" +str(Consensus), "-g" +str(Gaps), "-m" +str(minPeriod), "-x" +str(maxPeriod), "-a"+str(file), "-t1", "-T" +str(Coverage)]):
        print(file + " - XSTREAM run failed, skipping this file.")
        return False

    ## find html file which is output from XSTREAM, extract sequences, then run metapredict module against
    print(lineenter)

    html_files = [f for f in os.listdir(PATH) if f.endswith("_2.html")]
    if not html_files or not xstream_found_repeats(html_files):
        print("No tandem repeats detected by XSTREAM for " + str(file) + " - no candidates found, skipping.")
        print(lineenter)
        return False

    for f in html_files:
        ## extract sequence IDs, then retrieve sequences from input file and output temp file for iupred/filtering
        ## if XSTREAM filtering used, repeat region extract from _2.html file, filtered with input variables
        ## then sequence IDs used to extract sequences of filtered sequences from input file
        xstream_extract2(f, file, Aromatic, Electrostatic, metapredict_filter_value)
        ## remove the first XSTREAM run's files before regenerating the final report below
        for x in glob.glob("XSTREAM*"):
            os.remove(x)

        ## write out the final candidate set (already filtered above) + optional disorder plots
        metapredict_htp('Temp_xstream_filtered.fasta', directory, metapredict_plot)

        print("\nTidying up...")

        ## run XSTREAM once more, now on the final candidate set, so the saved HTML report
        ## reflects the actual final candidates (post-metapredict filtering) rather than an
        ## intermediate set. Skip it entirely if nothing survived the aromatic/electrostatic +
        ## disorder filters - XSTREAM crashes (ArrayIndexOutOfBoundsException) on an empty FASTA
        ## input, which would otherwise leave 0-byte report files behind with no indication of
        ## what happened.
        if not validate_fasta('candidate_sequences.fasta'):
            print("No candidates remained after aromatic/electrostatic and disorder filtering for " + str(file) + " - skipping final XSTREAM report.")
        elif not run_xstream(['java', '-Xmx1000m', '-Xms1000m', '-jar', XSTREAM_JAR, 'candidate_sequences.fasta', "-e"+str(Copy), "-i"+str(Word), "-I" +str(Consensus), "-g" +str(Gaps), "-m" +str(minPeriod), "-x" +str(maxPeriod), "-a"+str(file), "-t1", "-T" +str(Coverage)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT):
            print(file + " - final XSTREAM report run failed; candidate outputs may be incomplete.")
        else:
            ## build the interactive candidate report (see scripts/candidate_report.py), sourced
            ## from this final run's own "_2.html" - the candidates the user actually gets, not
            ## the earlier pre-metapredict-filtering set the first XSTREAM run above reported on
            final_html_files = [f for f in os.listdir(PATH) if f.endswith("_2.html")]
            if final_html_files and xstream_found_repeats(final_html_files):
                report_df, sequences, alignment_provider = parse_xstream_repeats(final_html_files[0], 'candidate_sequences.fasta')
                candidate_report_html = "{}_candidate_report.html".format(file)
                build_candidate_report(report_df, sequences, alignment_provider, candidate_report_html,
                                        score_meta={"label": "match quality", "explanation": XSTREAM_SCORE_EXPLANATION, "format": "{:.0f}%"})
                print("Candidate report written to " + candidate_report_html)

    return True

## module to characterize a FASTA file of known/reference target sequences (e.g. known pyrenoid
## linkers) and suggest FLIPPer search parameters from their observed pI, repeat structure and
## disorder. This is a standalone diagnostic action, not part of the proteome-scanning pipeline -
## it never filters anything, it just reports what the targets look like and suggests a starting
## point for --pi/--copy/--min-period/--max-period/--metapredict-filter-value.
def characterize(file):
    import os
    import glob
    import math
    import subprocess
    import pandas as pd
    import protfasta
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import metapredict as meta
    from FLIPPer_lib import lineenter

    print(lineenter)
    print("Characterizing target sequences in '" + file + "'\n")

    ## invalid_sequence_action='convert' strips/replaces characters like a trailing "*" stop
    ## codon, matching how the main pipeline's metapredict step (metapredict_htp) reads sequences
    records = protfasta.read_fasta(file, invalid_sequence_action='convert', return_list=True)
    if not records:
        print(file + " contains no valid sequences - aborting characterization.")
        return

    rows = []
    for seq_id, sequence in records:
        analysis = ProteinAnalysis(sequence)
        ## amino_acids_percent is 0-100 in current Biopython - divide by 100 to get the 0-1
        ## fraction that --serine/--alanine/--th-ratio thresholds are defined against (matching
        ## the same conversion in analysis_and_filtering)
        aa_percent = {aa: value / 100.0 for aa, value in analysis.amino_acids_percent.items()}
        helix = aa_percent['F'] + aa_percent['I'] + aa_percent['L'] + aa_percent['V'] + aa_percent['W'] + aa_percent['Y']
        turn = aa_percent['P'] + aa_percent['N'] + aa_percent['G'] + aa_percent['S']
        rows.append({
            'ID': seq_id,
            'Length': len(sequence),
            'pI': analysis.isoelectric_point(),
            'Serine %': aa_percent['S'],
            'Alanine %': aa_percent['A'],
            'TH Ratio': (turn / helix) if helix > 1e-9 else 10.0,
        })
    sequences = dict(records)

    ## --- XSTREAM repeat detection -------------------------------------------------------------

    ## runs XSTREAM once against the whole target file with the given extra flags, and returns
    ## {seq_id: [(start, end, period, copy_number), ...]} - every repeat region XSTREAM reported
    ## per sequence, not just the best one, so callers can decide what "still found" means
    def run_and_parse(extra_args):
        for x in glob.glob("XSTREAM*"):
            os.remove(x)
        xstream_args = ["java", "-Xmx1000m", "-Xms1000m", "-jar", XSTREAM_JAR, file] + extra_args + ["-a"+file, "-t1"]
        if not run_xstream(xstream_args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT):
            return {}
        html_files = glob.glob("XSTREAM_*_2.html")
        if not html_files:
            return {}
        single_id = records[0][0] if len(records) == 1 else None
        result = {}
        for row in parse_xstream_index(html_files[0], single_id):
            result.setdefault(row["seq_id"], []).append((row["begin"], row["end"], row["period"], row["copies"]))
        return result

    def overlaps_baseline(spans, baseline, min_frac=0.5):
        b_start, b_end = baseline
        baseline_len = b_end - b_start
        if baseline_len <= 0:
            return False
        for start, end, _, _ in spans:
            overlap = min(end, b_end) - max(start, b_start)
            if overlap / baseline_len >= min_frac:
                return True
        return False

    ## Baseline pass: deliberately loose copy/gaps/period-range/coverage - the goal here is to
    ## reliably find *a* repeat in each known-positive target, not to filter. Word/consensus match
    ## are held at the pipeline's own default (0.3625) rather than loosened further - looser than
    ## that lets XSTREAM lock onto a spurious short sub-harmonic of the real repeat (e.g. reporting
    ## a 31aa/11-copy unit instead of the real 65aa/5.3-copy one).
    baseline_report = run_and_parse(["-e2", "-i0.3625", "-I0.3625", "-g70", "-m10", "-x300", "-T0.2"])
    if not baseline_report:
        print("XSTREAM did not detect any tandem repeats in the target sequences - repeat statistics unavailable.")

    repeat_info = {}
    for seq_id, spans in baseline_report.items():
        start, end, period_val, copy_val = max(spans, key=lambda s: s[1] - s[0])
        repeat_info[seq_id] = {'Repeat Period (aa)': period_val, 'Repeat Copies': copy_val, 'span': (start, end)}

    ## disorder is reported over the detected repeat region itself where one was found (matching
    ## the main pipeline's now region-specific disorder filter), falling back to the whole protein
    ## for targets XSTREAM found no repeat in
    for row in rows:
        info = repeat_info.get(row['ID'], {})
        row['Repeat Period (aa)'] = info.get('Repeat Period (aa)')
        row['Repeat Copies'] = info.get('Repeat Copies')
        if 'span' in info:
            start, end = info['span']
            disorder_sequence = sequences[row['ID']][int(start) - 1:int(end)]
        else:
            disorder_sequence = sequences[row['ID']]
        row['Percent Disorder'] = meta.percent_disorder(disorder_sequence)

    ## --- word match / consensus match / gaps sensitivity sweep --------------------------------
    ## How permissive XSTREAM needs to be to still find each target's repeat varies a lot by
    ## target, and matters as much as period/copy for a real search. Rather than guess, sweep each
    ## parameter (holding the others at the loose baseline) and find the strictest setting that
    ## still detects a repeat overlapping the baseline span by at least half its length, for every
    ## target that had a baseline repeat. One XSTREAM run per grid point covers all targets at
    ## once, so this costs a fixed ~18 extra fast runs regardless of how many targets there are.
    baseline_spans = {seq_id: info['span'] for seq_id, info in repeat_info.items()}
    match_grid = [0.3625, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
    gaps_grid = [90, 80, 70, 60, 50, 40, 30, 20, 10, 0]
    max_match_tolerated = {}
    min_gaps_tolerated = {}
    if baseline_spans:
        for m in match_grid:
            report = run_and_parse(["-e2", "-i{:.4f}".format(m), "-I{:.4f}".format(m), "-g70", "-m10", "-x300", "-T0.2"])
            for seq_id, baseline in baseline_spans.items():
                if overlaps_baseline(report.get(seq_id, []), baseline):
                    max_match_tolerated[seq_id] = max(max_match_tolerated.get(seq_id, 0.0), m)
        for g in gaps_grid:
            report = run_and_parse(["-e2", "-i0.3625", "-I0.3625", "-g{}".format(g), "-m10", "-x300", "-T0.2"])
            for seq_id, baseline in baseline_spans.items():
                if overlaps_baseline(report.get(seq_id, []), baseline):
                    prior = min_gaps_tolerated.get(seq_id)
                    min_gaps_tolerated[seq_id] = g if prior is None else min(prior, g)

    df = pd.DataFrame(rows)

    ## clean up this exploratory XSTREAM run's files - they're not real pipeline output and would
    ## otherwise collide with a later real run against the same directory
    for x in glob.glob("XSTREAM*"):
        os.remove(x)

    no_extension = os.path.splitext(os.path.basename(file))[0]
    destination_folder = "{}_characterization".format(no_extension)
    os.makedirs(destination_folder, exist_ok=True)
    report_path = os.path.join(destination_folder, "{}_target_characterization.csv".format(no_extension))
    df.to_csv(report_path, index=None, sep=',')

    print(df.to_string(index=False))
    print("\nFull report written to " + report_path)

    missing = df[df['Repeat Period (aa)'].isna()]
    if len(missing):
        print("\nNo tandem repeat detected for: " + ", ".join(missing['ID'].tolist()) + " - excluded from repeat-based suggestions.")

    detected = df.dropna(subset=['Repeat Period (aa)', 'Repeat Copies'])

    def round_down_5(x):
        return 5 * math.floor(x / 5)

    def round_up_5(x):
        return 5 * math.ceil(x / 5)

    lines = [lineenter, "Suggested search parameters based on " + str(len(df)) + " target sequence(s):", ""]

    ## --pi is a minimum (pI >= threshold) OR maximum (pI <= threshold) filter depending on
    ## --pi-direction - the main pipeline defaults to "min" (basic/Arg-Lys-rich linkers like
    ## EPYC1/CsLinker). If the targets are consistently acidic instead, suggest "max" with a
    ## padded upper bound; if they straddle neutral pI, no single threshold/direction separates
    ## them from typical proteins, so flag that rather than suggest something misleadingly inert.
    pi_min, pi_max = df['pI'].min(), df['pI'].max()
    if pi_max < 7.0:
        pi_direction_suggest = "max"
        pi_suggest = round(min(14.0, pi_max + 1.0), 1)
        lines.append("\tpI: observed {:.2f} - {:.2f} (acidic)  ->  --pi {:.1f} --pi-direction max".format(pi_min, pi_max, pi_suggest))
    elif pi_min > 7.0:
        pi_direction_suggest = "min"
        pi_suggest = round(max(0.0, pi_min - 1.0), 1)
        lines.append("\tpI: observed {:.2f} - {:.2f} (basic)  ->  --pi {:.1f} --pi-direction min".format(pi_min, pi_max, pi_suggest))
    else:
        pi_direction_suggest = "min"
        pi_suggest = 0.0
        lines.append("\tpI: observed {:.2f} - {:.2f} - spans both acidic and basic, so no single --pi threshold/direction separates these targets from typical proteins. Suggesting --pi 0 (i.e. disabled); rely on repeat/disorder criteria instead.".format(pi_min, pi_max))

    ## --serine/--alanine/--th-ratio are minimum thresholds too - pad downward only. If the
    ## resulting suggestion is near zero, this trait isn't diagnostic for this target set (they
    ## simply aren't as Ser/Ala-rich or turn-dominated as EPYC1/CsLinker-type linkers), so say so
    ## rather than silently suggest a filter that will do essentially nothing.
    serine_suggest = round(max(0.0, df['Serine %'].min() - 0.02), 3)
    alanine_suggest = round(max(0.0, df['Alanine %'].min() - 0.01), 3)
    th_suggest = round(max(0.0, df['TH Ratio'].min() - 0.3), 2)
    lines.append("\tSerine content: observed {:.1%} - {:.1%}  ->  --serine {:.3f}{}".format(
        df['Serine %'].min(), df['Serine %'].max(), serine_suggest,
        "  (low - may not usefully separate these targets from typical proteins)" if serine_suggest <= 0.01 else ""))
    lines.append("\tAlanine content: observed {:.1%} - {:.1%}  ->  --alanine {:.3f}{}".format(
        df['Alanine %'].min(), df['Alanine %'].max(), alanine_suggest,
        "  (low - may not usefully separate these targets from typical proteins)" if alanine_suggest <= 0.002 else ""))
    lines.append("\tTurn/Helix ratio: observed {:.2f} - {:.2f}  ->  --th-ratio {:.2f}{}".format(
        df['TH Ratio'].min(), df['TH Ratio'].max(), th_suggest,
        "  (low - these targets aren't turn/coil-dominated the way EPYC1/CsLinker are; consider --th-ratio 0)" if th_suggest <= 0.2 else ""))

    disorder_suggest = max(0.0, round_down_5(df['Percent Disorder'].min() - 10.0))
    lines.append("\t% Disorder (of repeat region, or whole protein if none detected): observed {:.1f} - {:.1f}  ->  --metapredict-filter-value {:.0f}".format(
        df['Percent Disorder'].min(), df['Percent Disorder'].max(), disorder_suggest))

    example_flags = "--engine xstream --pi {:.1f} --pi-direction {} --th-ratio {:.2f} --serine {:.3f} --alanine {:.3f} --metapredict-filter-value {:.0f}".format(
        pi_suggest, pi_direction_suggest, th_suggest, serine_suggest, alanine_suggest, disorder_suggest)

    if len(detected):
        copy_suggest = max(2, math.floor(detected['Repeat Copies'].min() - 1))
        min_obs = detected['Repeat Period (aa)'].min()
        max_obs = detected['Repeat Period (aa)'].max()
        period_pad = max(10, round_down_5(0.15 * detected['Repeat Period (aa)'].median()))
        ## --min-period isn't a clean "only report periods >= this" cutoff - empirically, whether
        ## XSTREAM finds a given repeat at all as --min-period increases is NOT monotonic (e.g. a
        ## real period-56 repeat was found at -m20/25/30/40 but not at -m35 or -m42-45, with the
        ## *reported* period unchanged throughout). XSTREAM evidently detects periodicity via
        ## shorter internal seed comparisons and merges consecutive copies into the final repeat,
        ## and --min-period constrains that seed search too, not just the final answer - so a
        ## min-period comfortably below even the shortest observed period can still lose repeats
        ## unpredictably. The only reliable mitigation is to keep it generously low (close to
        ## FLIPPer's own default of 20) rather than tailored tightly to what the references show;
        ## only widen the window beyond the defaults (20-120) if a target needs it, never narrow
        ## inside them.
        DEFAULT_MIN_PERIOD, DEFAULT_MAX_PERIOD = 20.0, 120.0
        min_period_suggest = min(DEFAULT_MIN_PERIOD, max(10, round_down_5(min_obs - period_pad)))
        max_period_suggest = max(DEFAULT_MAX_PERIOD, round_up_5(max_obs + period_pad))
        lines.append("\tRepeat copies: observed {:.2f} - {:.2f}  ->  --copy {}".format(
            detected['Repeat Copies'].min(), detected['Repeat Copies'].max(), copy_suggest))
        lines.append("\tRepeat length: observed {:.0f} - {:.0f} aa  ->  --min-period {:.0f} --max-period {:.0f}{}".format(
            min_obs, max_obs, min_period_suggest, max_period_suggest,
            "  (kept at the pipeline defaults - narrower risks missing more divergent homologs)" if min_period_suggest == DEFAULT_MIN_PERIOD and max_period_suggest == DEFAULT_MAX_PERIOD else ""))
        example_flags += " --copy {} --min-period {:.0f} --max-period {:.0f}".format(copy_suggest, min_period_suggest, max_period_suggest)

        ## IMPORTANT: the sweep measures how well each target matches ITSELF - it's a ceiling on
        ## strictness, not a safe operating point. A real homolog you're searching for in a genome
        ## will be more divergent than the reference sequence is from itself, so tightening
        ## --word/--consensus/--gaps beyond the pipeline's own defaults just because the targets
        ## can tolerate it is likely to search too narrowly and return few or no real candidates.
        ## The sweep is genuinely useful in the other direction though: if even the default is
        ## stricter than a target can tolerate, that target needs loosening beyond default to be
        ## detected at all, and this catches that.
        DEFAULT_WORD, DEFAULT_CONSENSUS, DEFAULT_GAPS = 0.3625, 0.4, 55.0
        margin = 0.05
        if max_match_tolerated:
            ceiling = min(max_match_tolerated.values()) - margin
            limiting_id = min(max_match_tolerated, key=max_match_tolerated.get)
            word_suggest = round(min(DEFAULT_WORD, max(0.2, ceiling)), 4)
            consensus_suggest = round(min(DEFAULT_CONSENSUS, max(0.2, ceiling)), 4)
            if ceiling < DEFAULT_WORD:
                note = "  ({} needs looser matching than the pipeline default ({}) to be detected at all)".format(limiting_id, DEFAULT_WORD)
            else:
                note = "  (pipeline default already comfortably within every target's tolerance - not suggesting anything stricter than default, since that would only reflect self-similarity, not real sequence divergence)"
            lines.append("\tWord/consensus match tolerance (per target, self-similarity ceiling): " +
                          ", ".join("{}: up to {:.4f}".format(k, v) for k, v in sorted(max_match_tolerated.items())) +
                          "  ->  --word {:.4f} --consensus {:.4f}{}".format(word_suggest, consensus_suggest, note))
            example_flags += " --word {:.4f} --consensus {:.4f}".format(word_suggest, consensus_suggest)
        else:
            lines.append("\tWord/consensus match sweep found no tolerance data (targets may have failed at every grid point tested) - leaving --word/--consensus at their defaults.")

        if min_gaps_tolerated:
            floor = max(min_gaps_tolerated.values()) + 10
            limiting_id = max(min_gaps_tolerated, key=min_gaps_tolerated.get)
            gaps_suggest = max(DEFAULT_GAPS, floor)
            if floor > DEFAULT_GAPS:
                note = "  ({} needs more gaps allowed than the pipeline default ({:.0f}) to be detected at all)".format(limiting_id, DEFAULT_GAPS)
            else:
                note = "  (pipeline default already comfortably within every target's tolerance)"
            lines.append("\tGaps tolerance (per target, self-similarity floor): " +
                          ", ".join("{}: needs >= {}".format(k, v) for k, v in sorted(min_gaps_tolerated.items())) +
                          "  ->  --gaps {:.0f}{}".format(gaps_suggest, note))
            example_flags += " --gaps {:.0f}".format(gaps_suggest)
        else:
            lines.append("\tGaps sweep found no tolerance data (targets may have failed at every grid point tested) - leaving --gaps at its default.")

        lines.append("")
        lines.append("Example:")
        lines.append("\tpython3 FLIPPer.py --non-interactive " + example_flags)
        lines.append("")
        lines.append("These are a starting point, not a guarantee. --word/--consensus/--gaps above are never")
        lines.append("suggested stricter than FLIPPer's own defaults (0.3625/0.4/55), since a real search looks")
        lines.append("for divergent homologs, not copies of the targets themselves - so if a real search still")
        lines.append("returns few or no candidates, the more likely culprits are --pi/--th-ratio/--serine/")
        lines.append("--alanine (composition thresholds) or --min-period/--max-period/--coverage (repeat-shape")
        lines.append("constraints) being too strict for how divergent the real homologs are - try loosening")
        lines.append("those next.")
    else:
        lines.append("\tNo repeats detected in any target - cannot suggest --copy/--min-period/--max-period/--word/--consensus/--gaps; check the target sequences or loosen these manually.")
    lines.append(lineenter)

    report_text = "\n".join(lines)
    print(report_text)
    with open(os.path.join(destination_folder, "{}_suggested_parameters.txt".format(no_extension)), 'w') as f:
        f.write(report_text)

##module to output variables
def output_variables(file, pI, THRatio, Serine, Alanine, Copy, Word, Consensus, Gaps, minPeriod, maxPeriod, Coverage, Aromatic, Electrostatic, metapredict_filter_value, pI_direction="min"):
    import sys
    print("Outputting variables used.")
    tem = sys.stdout
    sys.stdout= m =open('{}_variables.txt'.format(file),'w')
    print("==========================================")
    print("Filtering variables:")
    print("\tpI Threshold: ", pI, "(direction: {}, i.e. keep pI {} threshold)".format(pI_direction, "<=" if str(pI_direction).lower() == "max" else ">="))
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

## engine-specific extra output artifacts finalize_output() (scripts/FLIPPer_lib.py) moves into
## the destination folder, beyond the common ones it already handles itself
def finalize_extra(destination_folder, path):
    import glob
    from FLIPPer_lib import move_files
    for z in glob.glob("XSTREAM*"):
        move_files(z, destination_folder, path)

## engine-specific temp files cleanup_temp_files() (scripts/FLIPPer_lib.py) removes, beyond the
## common one it already handles itself
def cleanup_extra():
    import os
    import glob
    for t in ["Temp_xstream_filtered.fasta", "Temp_Xstream_positives.fasta"]:
        if os.path.exists(t):
            os.remove(t)
    for x in glob.glob("XSTREAM*"):
        os.remove(x)
