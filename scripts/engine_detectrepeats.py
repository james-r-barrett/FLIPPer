## DECIPHER::DetectRepeats repeat-detection engine for FLIPPer.py - see engine_xstream.py for the
## other engine. Both modules expose the same small interface (NAME, EXTRA_PY_REQUIREMENTS,
## check_dependencies, process_file, characterize, output_variables, finalize_extra,
## cleanup_extra) so FLIPPer.py's shared driver can call either interchangeably via --engine.
import os

NAME = "detectrepeats"

## no extra Python packages beyond FLIPPer.py's shared requirements list
EXTRA_PY_REQUIREMENTS = []

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DETECT_REPEATS_R = os.path.join(SCRIPT_DIR, "detect_repeats.R")

## Default variables for this engine's own search parameters, overwritten by user input if y
## answered - used both directly (non-interactive mode with no override) and, via
## `from engine_detectrepeats import *`, as the interactive prompts' unmodified fallback
MinScore= "8"
MinCopies= "3"

## DetectRepeats' own maxPeriod argument is NOT a simple reporting ceiling - it constrains the
## periodicities the seed-and-extend search actually probes, and setting it close to (or below)
## a real repeat's period can make DetectRepeats miss that repeat entirely, even though the
## period is well under the value passed. (Confirmed empirically: CsLinker/EPYC1/SUPA1, true
## periods 45-71aa, were found with maxPeriod>=150 but NOT with maxPeriod=120.) So FLIPPer's
## user-facing --max-period is never passed into DetectRepeats itself - every run instead uses
## this constant (DetectRepeats' own documented default), and --max-period is enforced afterward
## as a post-filter on the reported Period column, the same way --min-period already has to be.
DETECT_REPEATS_SEARCH_MAX_PERIOD = 2000

## checked once, right after the engine is chosen, before asking the user anything else - every
## file's processing depends on Rscript and DECIPHER, and failing fast avoids wasting the user's
## time on prompts before a guaranteed failure
def check_dependencies():
    import subprocess
    try:
        rscript_ok = subprocess.call(["Rscript", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
    except FileNotFoundError:
        rscript_ok = False
    if not rscript_ok:
        print("Rscript could not be found on this system - DetectRepeats requires a working R installation.")
        return False
    try:
        decipher_ok = subprocess.call(
            ["Rscript", "-e", "suppressMessages(library(DECIPHER))"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        ) == 0
    except FileNotFoundError:
        decipher_ok = False
    if not decipher_ok:
        print("The R package DECIPHER could not be loaded - install it with BiocManager::install('DECIPHER').")
        return False
    return True

## minimum sequences per chunk when splitting work across parallel DetectRepeats processes -
## below this, the fixed ~1-2s cost of starting R and loading DECIPHER per process would
## dominate over the actual per-chunk work, so splitting further stops paying off.
DETECT_REPEATS_MIN_CHUNK_SIZE = 20

## module to run DECIPHER::DetectRepeats (via the detect_repeats.R wrapper), checking that it
## actually succeeded rather than assuming it did.
##
## DetectRepeats' own `processors` argument does not parallelize its tandem-repeat search
## (confirmed empirically: no speedup vs processors=1 on this DECIPHER version), so on a large
## candidate set this splits the input across up to os.cpu_count() separate Rscript processes -
## real, process-level parallelism - and merges their reports back into one CSV. Below
## DETECT_REPEATS_MIN_CHUNK_SIZE sequences it collapses to a single, unsplit Rscript call.
def run_detect_repeats(r_script, fasta_path, min_score, max_period, max_copies, out_csv, alignments_dir="", **kwargs):
    import subprocess
    import os
    import math
    import pandas as pd
    from Bio import SeqIO

    empty_columns = ['ID', 'Begin', 'End', 'Period', 'Copies', 'Score', 'RepeatIndex']

    def run_single(fasta, csv_out):
        args = ["Rscript", r_script, fasta, str(min_score), str(max_period), str(max_copies), csv_out, alignments_dir]
        try:
            result = subprocess.call(args, **kwargs)
        except FileNotFoundError as e:
            print("Could not run DetectRepeats: " + str(e))
            return False
        if result != 0:
            print("DetectRepeats exited with an error (code {}).".format(result))
            return False
        return True

    records = list(SeqIO.parse(fasta_path, "fasta-pearson"))
    if not records:
        pd.DataFrame(columns=empty_columns).to_csv(out_csv, index=False)
        return True

    n_chunks = max(1, min(os.cpu_count() or 1, math.ceil(len(records) / DETECT_REPEATS_MIN_CHUNK_SIZE)))
    if n_chunks == 1:
        return run_single(fasta_path, out_csv)

    chunk_size = math.ceil(len(records) / n_chunks)
    chunk_paths = []
    processes = []
    try:
        for i in range(n_chunks):
            chunk_records = records[i * chunk_size:(i + 1) * chunk_size]
            if not chunk_records:
                continue
            chunk_fasta = "{}.chunk{}.fasta".format(out_csv, i)
            chunk_csv = "{}.chunk{}.csv".format(out_csv, i)
            SeqIO.write(chunk_records, chunk_fasta, "fasta")
            chunk_paths.append((chunk_fasta, chunk_csv))
            ## all chunks share the same alignments_dir - safe to write concurrently since
            ## filenames are derived from each hit's own (globally unique) sequence ID
            args = ["Rscript", r_script, chunk_fasta, str(min_score), str(max_period), str(max_copies), chunk_csv, alignments_dir]
            try:
                processes.append(subprocess.Popen(args, **kwargs))
            except FileNotFoundError as e:
                print("Could not run DetectRepeats: " + str(e))
                return False

        if any(p.wait() != 0 for p in processes):
            print("DetectRepeats exited with an error in one or more parallel chunks.")
            return False

        frames = [pd.read_csv(csv_path) for _, csv_path in chunk_paths if os.path.exists(csv_path)]
        merged = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=empty_columns)
        merged.to_csv(out_csv, index=False)
        return True
    finally:
        for chunk_fasta, chunk_csv in chunk_paths:
            if os.path.exists(chunk_fasta):
                os.remove(chunk_fasta)
            if os.path.exists(chunk_csv):
                os.remove(chunk_csv)

## reads a FASTA file's sequences keyed BOTH by its raw header line (everything after '>', used
## verbatim by DECIPHER::DetectRepeats/Biostrings as each hit's ID) and by the short,
## whitespace-truncated first token (Biopython's record.id convention, used by FLIPPer's own
## Python-authored temp FASTAs like CandidateSequences_Temp.FASTA). DetectRepeats' report ID
## column follows whichever convention the FASTA it actually searched used - a header-only file
## has both conventions collapse to the same string, but a file with description text after the
## ID (e.g. metapredict_htp's candidate_sequences.fasta, built via protfasta which preserves the
## full header) only matches on the raw-header key. Keying on both means callers can look up a
## DetectRepeats report ID without needing to know in advance which convention applies, instead
## of the two silently mismatching and every row falling out of a Bio-.id-only lookup.
def _read_fasta_dual_keyed(path):
    sequences = {}
    header = None
    chunks = []

    def flush():
        if header is None:
            return
        seq = ''.join(chunks)
        sequences[header] = seq
        sequences[header.split()[0] if header.split() else header] = seq

    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                flush()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line.strip())
    flush()
    return sequences

## module to load a DetectRepeats report CSV (written by detect_repeats.R) as a dataframe,
## applying the search-window filters that DetectRepeats itself has no direct arguments for
## (minimum period, minimum copy number, minimum sequence coverage - DetectRepeats only exposes
## ceilings for period/copies via maxPeriod/maxCopies, and an overall significance threshold via
## minScore). Returns an empty dataframe (not an error) if the report has no rows or the file is
## missing/unreadable, since "DetectRepeats found nothing" is an expected outcome, not a failure.
def load_detect_repeats_report(csv_path, input_file, MinPeriod, MaxPeriod, MinCopies, Coverage):
    import os
    import pandas as pd
    columns = ['ID', 'Begin', 'End', 'Period', 'Copies', 'Score', 'RepeatIndex']
    if not os.path.exists(csv_path):
        return pd.DataFrame(columns=columns)
    df = pd.read_csv(csv_path)
    if df.empty:
        return pd.DataFrame(columns=columns)
    seq_lengths = {k: len(v) for k, v in _read_fasta_dual_keyed(input_file).items()}
    df['SeqLen'] = df['ID'].map(seq_lengths)
    df = df.dropna(subset=['SeqLen'])
    df['CoverageFraction'] = (df['End'] - df['Begin'] + 1) / df['SeqLen']
    filtered = df[
        (df['Period'] >= float(MinPeriod)) &
        (df['Period'] <= float(MaxPeriod)) &
        (df['Copies'] >= float(MinCopies)) &
        (df['CoverageFraction'] >= float(Coverage))
    ]
    return filtered

## module to extract candidate sequences from a DetectRepeats report, filtered to only include
## sequences with a repeat region that contains at least the requested number of aromatic and
## electrostatic residues AND is at least metapredict_filter_value% disordered. Since DetectRepeats
## reports plain Begin/End positions rather than a pre-computed consensus motif, the repeat region
## is sliced directly out of the original sequence rather than recovered by parsing a report. All
## three criteria are checked jointly against the same repeat instance - a sequence with several
## repeat regions qualifies if at least one of them satisfies all three together, not if different
## regions each satisfy a different criterion.
def detect_repeats_extract(csv_path, input_file, Aromatic, Electrostatic, MinPeriod, MaxPeriod, MinCopies, Coverage, metapredict_filter_value):
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import metapredict as meta
    from protfasta.utilities import convert_to_valid
    from FLIPPer_lib import lineenter

    filtered = load_detect_repeats_report(csv_path, input_file, MinPeriod, MaxPeriod, MinCopies, Coverage)
    number_before_filtering = filtered['ID'].nunique()

    sequences = {record.id: str(record.seq) for record in SeqIO.parse(input_file, "fasta-pearson")}

    positive_IDs = set()
    for _, row in filtered.iterrows():
        seq = sequences.get(row['ID'])
        if seq is None:
            continue
        ## Begin/End from DetectRepeats are 1-based inclusive positions
        region = seq[int(row['Begin']) - 1:int(row['End'])]
        Y = ProteinAnalysis(region)
        counts = Y.count_amino_acids()
        aromatic_count = counts['F'] + counts['W'] + counts['Y']
        electrostatic_count = counts['D'] + counts['E'] + counts['R'] + counts['K']
        ## real proteomes commonly carry ambiguous/ununsual residues (X, U, B, Z) - metapredict
        ## rejects anything but the 20 standard amino acids, so sanitize with the same conversion
        ## protfasta's invalid_sequence_action='convert' applies elsewhere in the pipeline
        region_clean = convert_to_valid(region.upper())
        if aromatic_count >= float(Aromatic) and electrostatic_count >= float(Electrostatic) \
                and region_clean and meta.percent_disorder(region_clean) > float(metapredict_filter_value):
            positive_IDs.add(row['ID'])

    number_after_filtering = 0
    with open("Temp_detectrepeats_filtered.fasta", "w") as f:
        for record in SeqIO.parse(input_file, "fasta-pearson"):
            if record.id in positive_IDs:
                SeqIO.write([record], f, "fasta")
                number_after_filtering += 1

    print(lineenter)
    print("Number of sequences with a qualifying repeat region (period/copies/coverage):", number_before_filtering)
    print("Number of sequences with a repeat region meeting "+str(Aromatic)+ ' aromatic (W/Y/F), ' +str(Electrostatic)+
          ' electrostatic (D/E/R/K) residues, and >'+str(metapredict_filter_value)+'% disorder (within that region):', number_after_filtering)
    print(lineenter)

## module to turn a raw DetectRepeats report (written by run_detect_repeats, one row per
## detected repeat region before period/copies/coverage filtering) into the final saved report -
## filtered down to the rows that actually meet the search window, in place of the unfiltered
## raw hits. Returns the number of qualifying repeat regions in the final report.
def finalize_detect_repeats_report(csv_path, input_file, MinPeriod, MaxPeriod, MinCopies, Coverage):
    filtered = load_detect_repeats_report(csv_path, input_file, MinPeriod, MaxPeriod, MinCopies, Coverage)
    filtered = filtered.rename(columns={'CoverageFraction': 'Coverage'}).drop(columns=['SeqLen'])
    filtered.to_csv(csv_path, index=None, sep=',')
    return len(filtered.index)

## same sanitisation detect_repeats.R applies to build alignment filenames from sequence IDs -
## must match exactly (character-for-character) so this can find the file R wrote for a given ID
def _detect_repeats_safe_id(seq_id):
    import re
    return re.sub(r'[^A-Za-z0-9_.\-]', '_', seq_id)

## DetectRepeats' Score is a log-odds statistical significance score: it comes from aligning the
## repeat copies against each other with a substitution matrix corrected for the sequence's own
## background amino acid composition, plus (with useEmpirical=TRUE, the default) added log-odds
## terms for how typical the repeat's copy number, unit length and composition are compared to a
## training set of known structural tandem repeats. Higher = more conserved / less likely to have
## arisen by chance - shown as a tooltip on each card's score pill rather than assumed knowledge.
DETECT_REPEATS_SCORE_EXPLANATION = (
    "Log-odds significance score (DECIPHER::DetectRepeats): repeat copies aligned against each "
    "other with a background-corrected substitution matrix, plus empirical log-odds terms for "
    "how typical this repeat's copy number/length/composition are versus known structural "
    "tandem repeats. Higher = more conserved / less likely by chance."
)

## the per-file "core": run DetectRepeats, extract/filter candidates, run metapredict, then
## rebuild the final DetectRepeats report + candidate_report.html against the actual final
## candidate set. Called from FLIPPer.py's shared per-file loop, right after
## analysis_and_filtering(). Returns True if processing reached the end (whether or not a final
## report could be built), False if it bailed out early (DetectRepeats run failed, or found no
## repeats at all) - the caller uses this the same way the old per-file loop used `continue`, to
## decide whether to print "Done!"/write the variables file.
def process_file(file, PATH, directory, metapredict_plot, metapredict_filter_value, Aromatic, Electrostatic,
                  MinScore, MinCopies, minPeriod, maxPeriod, Coverage):
    import os
    import subprocess
    import pandas as pd
    from FLIPPer_lib import lineenter, metapredict_htp, build_candidate_report

    ## run DetectRepeats using filtered sequences. max_period/max_copies passed to R are
    ## always the generous DETECT_REPEATS_SEARCH_MAX_PERIOD/1000 - NOT the user's
    ## --min-period/--max-period/--min-copies window, which is enforced afterward in Python
    ## (see DETECT_REPEATS_SEARCH_MAX_PERIOD's definition for why max_period in particular
    ## can't just be the user's ceiling).
    candidates_csv = "DetectRepeats_candidates.csv"
    if not run_detect_repeats(DETECT_REPEATS_R, "CandidateSequences_Temp.FASTA", MinScore, DETECT_REPEATS_SEARCH_MAX_PERIOD, 1000, candidates_csv):
        print(file + " - DetectRepeats run failed, skipping this file.")
        return False

    print(lineenter)

    raw_report = pd.read_csv(candidates_csv) if os.path.exists(candidates_csv) else pd.DataFrame()
    if raw_report.empty:
        print("No tandem repeats detected by DetectRepeats for " + str(file) + " - no candidates found, skipping.")
        print(lineenter)
        return False

    ## extract sequence IDs whose repeat region meets the min-period/min-copies/coverage
    ## window, the aromatic/electrostatic composition filter, and the metapredict disorder
    ## filter - all three checked against that same repeat region
    detect_repeats_extract(candidates_csv, file, Aromatic, Electrostatic, minPeriod, maxPeriod, MinCopies, Coverage, metapredict_filter_value)
    os.remove(candidates_csv)

    ## write out the final candidate set (already filtered above) + optional disorder plots
    metapredict_htp('Temp_detectrepeats_filtered.fasta', directory, metapredict_plot)

    print("\nTidying up...")

    ## run DetectRepeats once more, now on the final candidate set, so the saved report
    ## reflects the actual final candidates (post-metapredict filtering) rather than an
    ## intermediate set. alignments_dir makes this run also extract and align each hit's
    ## individual repeat copies (used below to build the candidate report) - only done here,
    ## on the small final candidate set, not on the first (much larger) pre-filter pass.
    final_report_csv = "{}_detected_repeats.csv".format(file)
    alignments_dir = "DetectRepeats_alignments"
    if not run_detect_repeats(DETECT_REPEATS_R, 'candidate_sequences.fasta', MinScore, DETECT_REPEATS_SEARCH_MAX_PERIOD, 1000, final_report_csv,
                               alignments_dir=alignments_dir, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT):
        print(file + " - final DetectRepeats report run failed; candidate outputs may be incomplete.")
    elif os.path.exists(final_report_csv):
        finalize_detect_repeats_report(final_report_csv, 'candidate_sequences.fasta', minPeriod, maxPeriod, MinCopies, Coverage)

        from Bio import SeqIO
        report_df = pd.read_csv(final_report_csv) if os.path.exists(final_report_csv) else pd.DataFrame()
        sequences = _read_fasta_dual_keyed('candidate_sequences.fasta')

        def alignment_provider(seq_id, repeat_index):
            align_path = os.path.join(alignments_dir, "{}__{}.fasta".format(_detect_repeats_safe_id(seq_id), repeat_index))
            if not os.path.exists(align_path):
                return None
            return [(r.id, str(r.seq)) for r in SeqIO.parse(align_path, "fasta")]

        candidate_report_html = "{}_candidate_report.html".format(file)
        build_candidate_report(report_df, sequences, alignment_provider, candidate_report_html,
                                score_meta={"label": "score", "explanation": DETECT_REPEATS_SCORE_EXPLANATION, "format": "{:.1f}"})
        print("Candidate report written to " + candidate_report_html)

    return True

## module to characterize a FASTA file of known/reference target sequences (e.g. known pyrenoid
## linkers) and suggest FLIPPer search parameters from their observed pI, repeat structure and
## disorder. This is a standalone diagnostic action, not part of the proteome-scanning pipeline -
## it never filters anything, it just reports what the targets look like and suggests a starting
## point for --pi/--serine/--alanine/--th-ratio/--min-copies/--min-period/--max-period/--coverage/
## --min-score/--aromatic/--electrostatic/--metapredict-filter-value.
def characterize(file):
    import os
    import glob
    import math
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

    ## run DetectRepeats with a deliberately loose minScore - the goal here is to reliably find
    ## *a* repeat in each known-positive target, not to filter, so it's far more permissive than
    ## the main pipeline's default. max_period is always DETECT_REPEATS_SEARCH_MAX_PERIOD (see its
    ## definition above) - there is no "loose" vs "strict" max_period, only a value the search
    ## itself can actually work with.
    report_csv = "DetectRepeats_characterize_report.csv"
    if os.path.exists(report_csv):
        os.remove(report_csv)
    repeat_info = {}
    if not run_detect_repeats(DETECT_REPEATS_R, file, min_score=4, max_period=DETECT_REPEATS_SEARCH_MAX_PERIOD, max_copies=1000, out_csv=report_csv):
        print("DetectRepeats failed to run against " + file + " - repeat statistics unavailable; pI/disorder still reported.")
    elif not os.path.exists(report_csv):
        print("DetectRepeats produced no report for " + file + " - repeat statistics unavailable.")
    else:
        report_df = pd.read_csv(report_csv)
        if report_df.empty:
            print("DetectRepeats did not detect any tandem repeats in the target sequences - repeat statistics unavailable.")
        else:
            ## a sequence can have more than one repeat region reported; keep the one with the
            ## largest position span (the main repeat), not just whichever row came last
            report_df['Span'] = report_df['End'] - report_df['Begin']
            for seq_id, group in report_df.groupby('ID'):
                best = group.loc[group['Span'].idxmax()]
                repeat_info[seq_id] = {'Repeat Period (aa)': best['Period'], 'Repeat Copies': best['Copies'],
                                       'Repeat Score': best['Score'],
                                       'Begin': int(best['Begin']), 'End': int(best['End'])}

    if os.path.exists(report_csv):
        os.remove(report_csv)

    ## disorder/aromatic/electrostatic are reported over the detected repeat region itself where
    ## one was found (matching the main pipeline's region-specific filters in
    ## detect_repeats_extract), falling back to the whole protein for targets DetectRepeats found
    ## no repeat in - Serine/Alanine/TH Ratio, by contrast, are whole-protein composition figures
    ## the same way the main pipeline's own pre-repeat-detection filter (analysis_and_filtering)
    ## computes them, so they're never region-restricted
    rows = []
    for seq_id, sequence in records:
        info = repeat_info.get(seq_id, {})
        region = sequence[info['Begin'] - 1:info['End']] if 'Begin' in info else sequence
        coverage = (info['End'] - info['Begin'] + 1) / len(sequence) if 'Begin' in info else None
        analysis = ProteinAnalysis(sequence)
        ## amino_acids_percent is 0-100 in current Biopython - divide by 100 to get the 0-1
        ## fraction --serine/--alanine/--th-ratio are defined against, matching the same
        ## conversion analysis_and_filtering applies to the real pipeline's own pre-filter
        aa_percent = {aa: value / 100.0 for aa, value in analysis.amino_acids_percent.items()}
        helix = aa_percent['F'] + aa_percent['I'] + aa_percent['L'] + aa_percent['V'] + aa_percent['W'] + aa_percent['Y']
        turn = aa_percent['P'] + aa_percent['N'] + aa_percent['G'] + aa_percent['S']
        region_counts = ProteinAnalysis(region).count_amino_acids()
        rows.append({
            'ID': seq_id,
            'Length': len(sequence),
            'pI': analysis.isoelectric_point(),
            'Serine %': aa_percent['S'],
            'Alanine %': aa_percent['A'],
            'TH Ratio': (turn / helix) if helix > 1e-9 else 10.0,
            'Percent Disorder': meta.percent_disorder(region),
            'Repeat Period (aa)': info.get('Repeat Period (aa)'),
            'Repeat Copies': info.get('Repeat Copies'),
            'Repeat Score': info.get('Repeat Score'),
            'Repeat Coverage': coverage,
            'Repeat Aromatic': region_counts['F'] + region_counts['W'] + region_counts['Y'],
            'Repeat Electrostatic': region_counts['D'] + region_counts['E'] + region_counts['R'] + region_counts['K'],
        })

    df = pd.DataFrame(rows)

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

    def round_down_frac(x, step=0.05):
        return step * math.floor(x / step)

    lines = [lineenter, "Suggested search parameters based on " + str(len(df)) + " target sequence(s):", ""]

    ## --pi/--metapredict-filter-value/--serine/--alanine/--th-ratio are minimum thresholds in the
    ## main pipeline, so they're padded downward only - padding both ways would suggest excluding
    ## the targets themselves. --serine/--alanine/--th-ratio in particular are enforced by
    ## analysis_and_filtering() BEFORE repeat detection even runs - a target failing one of these
    ## never reaches DetectRepeats at all, so they're computed here unconditionally (whole-protein
    ## composition, not repeat-region-dependent) rather than skipped the way repeat-derived
    ## suggestions below are when nothing was detected.
    pi_suggest = round(max(0.0, df['pI'].min() - 1.0), 1)
    serine_suggest = round(max(0.0, df['Serine %'].min() - 0.02), 3)
    alanine_suggest = round(max(0.0, df['Alanine %'].min() - 0.01), 3)
    th_suggest = round(max(0.0, df['TH Ratio'].min() - 0.3), 2)
    disorder_suggest = max(0.0, round_down_5(df['Percent Disorder'].min() - 10.0))
    lines.append("\tpI: observed {:.2f} - {:.2f}  ->  --pi {:.1f}".format(df['pI'].min(), df['pI'].max(), pi_suggest))
    lines.append("\tSerine content: observed {:.1%} - {:.1%}  ->  --serine {:.3f}{}".format(
        df['Serine %'].min(), df['Serine %'].max(), serine_suggest,
        "  (low - may not usefully separate these targets from typical proteins)" if serine_suggest <= 0.01 else ""))
    lines.append("\tAlanine content: observed {:.1%} - {:.1%}  ->  --alanine {:.3f}{}".format(
        df['Alanine %'].min(), df['Alanine %'].max(), alanine_suggest,
        "  (low - may not usefully separate these targets from typical proteins)" if alanine_suggest <= 0.002 else ""))
    lines.append("\tTurn/Helix ratio: observed {:.2f} - {:.2f}  ->  --th-ratio {:.2f}{}".format(
        df['TH Ratio'].min(), df['TH Ratio'].max(), th_suggest,
        "  (low - these targets aren't turn/coil-dominated the way EPYC1/CsLinker are; consider --th-ratio 0)" if th_suggest <= 0.2 else ""))
    lines.append("\t% Disorder (of repeat region, or whole protein if none detected): observed {:.1f} - {:.1f}  ->  --metapredict-filter-value {:.0f}".format(
        df['Percent Disorder'].min(), df['Percent Disorder'].max(), disorder_suggest))

    example_flags = "--pi {:.1f} --serine {:.3f} --alanine {:.3f} --th-ratio {:.2f} --metapredict-filter-value {:.0f}".format(
        pi_suggest, serine_suggest, alanine_suggest, th_suggest, disorder_suggest)

    if len(detected):
        copy_suggest = max(2, math.floor(detected['Repeat Copies'].min() - 1))
        period_pad = max(10, round_down_5(0.15 * detected['Repeat Period (aa)'].median()))
        min_period_suggest = max(10, round_down_5(detected['Repeat Period (aa)'].min() - period_pad))
        max_period_suggest = round_up_5(detected['Repeat Period (aa)'].max() + period_pad)
        ## --coverage is also a minimum threshold (main pipeline keeps hits with
        ## CoverageFraction >= --coverage) - pad downward like pi/copies/disorder above. Without
        ## this, targets whose repeat region only spans part of the sequence (e.g. flanking
        ## transit/signal peptide, or linker/terminal residues outside the repeat) fall below the
        ## pipeline's 0.75 default and get silently dropped even though characterize found their
        ## repeat just fine.
        coverage_suggest = max(0.0, round(round_down_frac(detected['Repeat Coverage'].min() - 0.05), 2))
        ## --min-score is NOT a post-hoc filter like the ones above - it's passed straight into
        ## DECIPHER::DetectRepeats() as a hard search-time cutoff (see detect_repeats.R), the same
        ## way DETECT_REPEATS_SEARCH_MAX_PERIOD's own docstring warns maxPeriod is. A repeat scoring
        ## below --min-score is never found at all, not merely filtered out afterward - it won't
        ## even appear in the raw candidates CSV, unlike a coverage/period/copies miss. Since the
        ## consequence of guessing too high is silent and much harder to diagnose than guessing too
        ## low, this uses a larger relative pad (25%, floor 2) than the other minimum thresholds
        ## here rather than a small fixed offset.
        score_pad = max(2.0, 0.25 * detected['Repeat Score'].min())
        min_score_suggest = max(0, math.floor(detected['Repeat Score'].min() - score_pad))
        ## --aromatic/--electrostatic are also minimum thresholds, checked against the same
        ## detected repeat region (detect_repeats_extract) - pad downward the same way, with a
        ## floor of 1 residue so a target that just barely clears the observed minimum doesn't get
        ## a suggestion of 0 padding down to nothing.
        aromatic_pad = max(1, round(0.2 * detected['Repeat Aromatic'].min()))
        aromatic_suggest = max(0, math.floor(detected['Repeat Aromatic'].min() - aromatic_pad))
        electrostatic_pad = max(1, round(0.2 * detected['Repeat Electrostatic'].min()))
        electrostatic_suggest = max(0, math.floor(detected['Repeat Electrostatic'].min() - electrostatic_pad))
        lines.append("\tRepeat copies: observed {:.2f} - {:.2f}  ->  --min-copies {}".format(
            detected['Repeat Copies'].min(), detected['Repeat Copies'].max(), copy_suggest))
        lines.append("\tRepeat length: observed {:.0f} - {:.0f} aa  ->  --min-period {} --max-period {}".format(
            detected['Repeat Period (aa)'].min(), detected['Repeat Period (aa)'].max(), min_period_suggest, max_period_suggest))
        lines.append("\tRepeat region coverage of full sequence: observed {:.2f} - {:.2f}  ->  --coverage {:.2f}".format(
            detected['Repeat Coverage'].min(), detected['Repeat Coverage'].max(), coverage_suggest))
        lines.append("\tRepeat significance score (DetectRepeats, this run's own loose min-score=4 probe): observed {:.1f} - {:.1f}  ->  --min-score {}".format(
            detected['Repeat Score'].min(), detected['Repeat Score'].max(), min_score_suggest))
        lines.append("\tAromatic residues (W/Y/F) in repeat region: observed {:.0f} - {:.0f}  ->  --aromatic {}".format(
            detected['Repeat Aromatic'].min(), detected['Repeat Aromatic'].max(), aromatic_suggest))
        lines.append("\tElectrostatic residues (D/E/R/K) in repeat region: observed {:.0f} - {:.0f}  ->  --electrostatic {}".format(
            detected['Repeat Electrostatic'].min(), detected['Repeat Electrostatic'].max(), electrostatic_suggest))
        example_flags += " --min-copies {} --min-period {} --max-period {} --coverage {:.2f} --min-score {} --aromatic {} --electrostatic {}".format(
            copy_suggest, min_period_suggest, max_period_suggest, coverage_suggest, min_score_suggest, aromatic_suggest, electrostatic_suggest)
        lines.append("")
        lines.append("Example:")
        lines.append("\tpython3 FLIPPer.py --non-interactive --engine detectrepeats " + example_flags)
    else:
        lines.append("\tNo repeats detected in any target - cannot suggest --min-copies/--min-period/--max-period/--coverage/--min-score/--aromatic/--electrostatic; check the target sequences or loosen these manually.")
    lines.append(lineenter)

    report_text = "\n".join(lines)
    print(report_text)
    with open(os.path.join(destination_folder, "{}_suggested_parameters.txt".format(no_extension)), 'w') as f:
        f.write(report_text)

##module to output variables
def output_variables(file, pI, THRatio, Serine, Alanine, MinScore, MinCopies, minPeriod, maxPeriod, Coverage, Aromatic, Electrostatic, metapredict_filter_value, pI_direction="min"):
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
    print("DetectRepeats variables:")
    print("\tMinimum score: ", MinScore)
    print("\tMinimum copy number: ", MinCopies)
    print("\tMinimum repeat period: ", minPeriod)
    print("\tMaximum repeat period: ", maxPeriod)
    print("\tSequence coverage: ", Coverage)
    print("==========================================")
    print("Post-DetectRepeats Filtering")
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
    for z in glob.glob("*_detected_repeats.csv"):
        move_files(z, destination_folder, path)

## engine-specific temp files cleanup_temp_files() (scripts/FLIPPer_lib.py) removes, beyond the
## common one it already handles itself
def cleanup_extra():
    import os
    import glob
    for t in ["Temp_detectrepeats_filtered.fasta", "DetectRepeats_candidates.csv", "DetectRepeats_characterize_report.csv"]:
        if os.path.exists(t):
            os.remove(t)
    ## belt-and-braces sweep for parallel-chunk temp files from run_detect_repeats() - normally
    ## removed by its own finally block, this only matters if the process was killed mid-run
    for x in glob.glob("*.chunk*.fasta") + glob.glob("*.chunk*.csv"):
        os.remove(x)
    ## per-hit repeat-unit alignments are scratch input to the candidate report - once it's been
    ## built (or the run failed before getting that far), these aren't needed any more
    if os.path.exists("DetectRepeats_alignments"):
        import shutil
        shutil.rmtree("DetectRepeats_alignments")
