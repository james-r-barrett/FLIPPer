<p align="center"><img src="icon.svg" width="160" alt="FLIPPer icon"></p>

# 🔬 Fast Linker Identification Pipeline for Pyrenoids (FLIPPer)

FLIPPer runs on **FASTA-formatted protein sequences** to identify candidate linker regions associated with pyrenoids.

By default FLIPPer scans its **input directory** (the current directory, or wherever `--input-dir` points) for valid FASTA files and processes each one using the defined pipeline settings.

> **Reproducing published results:** the FLIPPer version used for our published linker searches is preserved as the [`v1.1`](../../releases/tag/v1.1) tag/release. Everything below documents the current `--engine`-based pipeline (v3.0), which changed default parameters and internals since then - check out `v1.1` instead if you need to exactly reproduce those results.

FLIPPer's tandem-repeat detection step is pluggable - pick the engine with `--engine`:

- **`xstream`** (default) - [XSTREAM](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-382), a Java tool bundled with FLIPPer. Requires a working **Java** installation.
- **`detectrepeats`** - [DetectRepeats](https://doi.org/10.1093/nar/gkaf866), from the R/Bioconductor package **DECIPHER**. Requires a working **R** installation with DECIPHER installed.

Both engines are filtered/scored through the same downstream pipeline (pI, composition, disorder) and produce an identical-looking interactive candidate report - only the repeat-detection step itself, and its own set of search parameters, differs.

---

## 📥 Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/james-r-barrett/FLIPPer.git
   cd FLIPPer
   ```

2. **Install the Python dependencies** (requires Python ≥ 3.10):

   ```bash
   pip install -r requirements.txt
   ```

   Alternatively, install them individually:

   ```bash
   pip install biopython matplotlib beautifulsoup4 pandas metapredict cython protfasta
   ```

3. **Install the dependencies for whichever repeat-detection engine(s) you plan to use:**

   - **`--engine xstream`** (default): requires a working **Java** installation. XSTREAM itself (`scripts/xstream.jar`) is bundled with FLIPPer, so there's nothing else to install.
   - **`--engine detectrepeats`**: requires **R** with the Bioconductor package **DECIPHER** installed:

     ```r
     if (!requireNamespace("BiocManager", quietly = TRUE))
         install.packages("BiocManager")
     BiocManager::install("DECIPHER")
     ```

That's it - FLIPPer is installed. You can run it in place from the cloned directory, or point `--input-dir` at data anywhere else (see [Usage](#-usage) below) without needing to copy `FLIPPer.py` alongside it.

---

## 🚀 Usage

```bash
python3 FLIPPer.py
```

Running this command from the terminal or command prompt processes every FASTA file in the current directory using **XSTREAM with the documented default parameters** below - no prompts, nothing to answer. This is the mode most users want most of the time.

FLIPPer can be installed once and run against FASTA files in any directory - it no longer needs to be copied alongside your data. Point it at your data with `--input-dir`:

```bash
python3 /path/to/FLIPPer/FLIPPer.py --input-dir /path/to/your/fasta_files
```

Override individual parameters with flags, e.g.:

```bash
python3 FLIPPer.py --input-dir /path/to/your/fasta_files --engine detectrepeats --pi 7 --coverage 0.4
```

Run `python3 FLIPPer.py --help` for the full list of flags and their defaults; each engine-specific flag's help text notes which `--engine` it applies to.

XSTREAM engine, every parameter spelled out:

```bash
python3 FLIPPer.py --engine xstream --input-dir /path/to/your/fasta_files \
  --pi 8 --pi-direction min --th-ratio 1 --serine 0.05 --alanine 0.01 \
  --copy 3 --word 0.3625 --consensus 0.4 --gaps 55 --min-period 20 --max-period 120 --coverage 0.75 \
  --aromatic 1 --electrostatic 2 \
  --metapredict-filter-value 50 --plots
```

DetectRepeats engine, every parameter spelled out:

```bash
python3 FLIPPer.py --engine detectrepeats --input-dir /path/to/your/fasta_files \
  --pi 8 --pi-direction min --th-ratio 1 --serine 0.05 --alanine 0.01 \
  --min-score 8 --min-copies 3 --min-period 20 --max-period 120 --coverage 0.75 \
  --aromatic 1 --electrostatic 2 \
  --metapredict-filter-value 50 --plots
```

---

### 🧭 Interactive mode

To be walked through every parameter instead, pass `--interactive`. Each prompt shows its current value (the flag you passed, or its documented default) in brackets - press Enter to keep it, or type a replacement:

```bash
python3 FLIPPer.py --input-dir /path/to/your/fasta_files --interactive
```

```
Which repeat-detection engine? (xstream/detectrepeats) [xstream]:
pI threshold [8.0]:
pI direction - keep pI >= threshold ('min', ...) or pI <= threshold ('max', ...) (min/max) [min]:
...
```

You can combine this with flags to change what the bracketed defaults are - e.g. `--interactive --pi 7` starts the pI prompt at `[7.0]` instead of `[8.0]`.

**The two engines' search parameters aren't the same set.** XSTREAM exposes its own search knobs directly (`--copy`/`--word`/`--consensus`/`--gaps`). DetectRepeats has no minimum copy number, word-match, consensus-match, or gap-count arguments - only an overall significance threshold (`--min-score`) plus internal ceilings on period/copy number - so FLIPPer enforces `--min-period`/`--max-period`/`--min-copies`/`--coverage` itself as post-filters on DetectRepeats' reported hits, keeping the same search-window behaviour XSTREAM gets from its own native flags.

---

### 🎯 Characterizing target/reference sequences

If you have a FASTA file of known linkers (e.g. EPYC1, CsLinker) and want a starting point for your search parameters, run:

```bash
python3 FLIPPer.py --engine xstream --characterize targets.fasta
```

`example_linkers.fasta`, bundled with FLIPPer, is one such file (EPYC1, CsLinker, SUPA1) - try the command above on it directly:

```bash
python3 FLIPPer.py --engine xstream --characterize example_linkers.fasta
```

For each target this reports pI, Serine/Alanine content, Turn/Helix ratio, repeat length/copy number, and % disorder (metapredict), then suggests search parameters padded around what was observed - for `--engine xstream`: `--pi`/`--pi-direction`, `--serine`, `--alanine`, `--th-ratio`, `--copy`, `--min-period`/`--max-period`, `--metapredict-filter-value`, and `--word`/`--consensus`/`--gaps`; for `--engine detectrepeats`: `--pi`, `--min-copies`, `--min-period`/`--max-period`, and `--metapredict-filter-value`. It runs standalone - it doesn't touch `--input-dir` or run the main pipeline - and writes its report to `<targets>_characterization/` alongside the input file.

A few things worth knowing about how the XSTREAM engine derives suggestions:

- **`--pi`/`--pi-direction`**: FLIPPer's pI filter can keep proteins *at least as basic as* a threshold (`min`, the default - suited to Arg/Lys-rich linkers like EPYC1/CsLinker) or *at least as acidic as* a threshold (`max`). The characterizer looks at the targets' pI and suggests whichever direction fits; if the targets span both acidic and basic pI, no single threshold/direction separates them from typical proteins, and it says so rather than suggest something misleadingly inert.
- **`--serine`/`--alanine`/`--th-ratio`**: if a suggested value comes out near zero, that trait isn't distinguishing for this target set (they simply aren't as Ser/Ala-rich or turn-dominated as EPYC1/CsLinker), and the report flags this so you know to loosen or drop that filter rather than trust a threshold that will pass almost everything.
- **`--word`/`--consensus`/`--gaps`**: these strongly affect whether XSTREAM detects a given repeat at all. The characterizer empirically sweeps a grid of values per target (not just a single fixed run) to find the strictest match/consensus threshold and fewest gaps each target still tolerates. That tolerance reflects how well a target matches *itself*, though, not how divergent a real homolog elsewhere in a proteome will be - so the suggestion is never tightened beyond FLIPPer's own defaults (0.3625/0.4/55), only loosened below them if a target needs it. This adds ~15-30 seconds of extra (fast) XSTREAM runs per characterization.
- **`--min-period`/`--max-period`**: `--min-period` is not a clean "only report periods >= this" cutoff - empirically, whether XSTREAM finds a real repeat at all as `--min-period` increases is *not monotonic* (a real period-56 repeat was found with `--min-period 20/25/30/40` but not `35` or `42-45`, with the period XSTREAM reported unchanged throughout). XSTREAM evidently detects periodicity via shorter internal seed comparisons and merges consecutive copies into the final repeat, and `--min-period` constrains that seed search too, not just the reported answer - so a value comfortably below even the shortest reference period can still lose real repeats unpredictably. The characterizer therefore never suggests narrowing this window below FLIPPer's own defaults (20-120) - only widening it if a target's observed period falls outside them.

**These suggestions are a starting point, not a guarantee.** If a real search using them returns few or no candidates, the likely culprits are the composition thresholds (`--pi`/`--th-ratio`/`--serine`/`--alanine`) or the repeat-shape constraints (`--min-period`/`--max-period`/`--coverage`) being too tightly scoped to the exact reference sequences - try loosening those, or falling back to the plain pipeline defaults as a baseline, before assuming there's nothing there.

---

### 🔁 Re-filtering an existing `--engine detectrepeats` run

Most of DetectRepeats' search-window filters (`--min-period`/`--max-period`/`--min-copies`/`--coverage`/`--aromatic`/`--electrostatic`/`--metapredict-filter-value`) are applied *after* detection, against a raw report of every candidate repeat DetectRepeats found. FLIPPer saves that raw report (`<file>_raw_repeats.csv`, plus a small sidecar) alongside its normal outputs, so if your first run's thresholds were too strict (or too loose) you can try different ones without re-running detection on the whole input:

```bash
python3 FLIPPer.py --engine detectrepeats --refilter path/to/file.fasta_FLIPPer_outputs \
  --coverage 0.4 --min-copies 4
```

This re-applies the given thresholds to the saved raw report, re-runs the (usually much smaller) surviving candidate set through DetectRepeats once more to regenerate alignments, and rebuilds the candidate report - writing everything into `<OUTPUT_DIR>/refiltered/` without touching the original run's own outputs. It runs standalone, like `--characterize` - it doesn't touch `--input-dir` or run the main pipeline.

**`--min-score` can't be changed this way.** Unlike the other filters, it's a search-time cutoff DetectRepeats itself applies - a repeat scoring below it is never written to the raw report in the first place, so there's nothing to re-filter. Changing it requires a full re-run.

---

## 📂 Outputs

For each analysed FASTA file, FLIPPer generates the following outputs:

### 1️⃣ Metapredict Plots
- PDF graphs containing:
  - Metapredict disorder scores  
  - Predicted AlphaFold2 pLDDT values  

---

### 2️⃣ Candidate Report
- A single self-contained `*_candidate_report.html` - one collapsible card per final candidate, combining its sequence/properties, an interactive disorder + predicted pLDDT chart with each repeat region plotted against it, and its repeat copies aligned against a consensus
- Identical in format between both engines - only the score pill's meaning differs (XSTREAM: match quality, derived from ConsensusError; DetectRepeats: log-odds significance score), explained in the pill's tooltip

---

### 3️⃣ Engine-specific repeat detection output
- **`--engine xstream`**: XSTREAM's own tandem repeat detection results, viewable as a linked HTML file (`*_out_2.html`)
- **`--engine detectrepeats`**: a CSV report, `<file>_detected_repeats.csv`, one row per detected repeat region (`ID`, `Begin`, `End`, `Period`, `Copies`, `Score`, `Coverage`) - plus `<file>_raw_repeats.csv` (every candidate repeat found, before the period/copies/coverage/aromatic/electrostatic/disorder filters) and a small `.meta.json` sidecar, used by [`--refilter`](#-re-filtering-an-existing---engine-detectrepeats-run) to re-apply different thresholds without re-running detection

---

### 4️⃣ Candidate Sequences
- FASTA file containing all candidate linker sequences  
- CSV file containing the same sequences  

---

### 5️⃣ Variables File
- Text file containing the user-defined variables used during the pipeline run  

---

### 6️⃣ Sequence Analysis
- CSV files containing sequence analysis results generated during Step 1 of sequence filtering  

---

## 📦 Requirements

- Python ≥ 3.10
- **`--engine xstream`**: Java (bundled `xstream.jar` requires a working Java installation)
- **`--engine detectrepeats`**: R with the Bioconductor package `DECIPHER` installed
- biopython
- matplotlib
- pandas
- metapredict
- cython
- protfasta
- beautifulsoup4 (only needed for `--engine xstream`)

See [Installation](#-installation) above for how to install all of these from a fresh checkout.

---

## 🧬 Third-Party Software & Citations

### XSTREAM

This pipeline utilises **XSTREAM**, a tandem repeat detection algorithm developed by Aaron Newman and James Cooper.

If you use FLIPPer with `--engine xstream`, please cite:

Newman, A. M. & Cooper, J. B. (2007).  
*XSTREAM: a practical algorithm for identification and architecture modeling of tandem repeats in protein sequences.*  
BMC Bioinformatics, 8, 382.  
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-382

---

### DetectRepeats / DECIPHER

This pipeline optionally utilises **DetectRepeats**, a tandem repeat detection method in the R/Bioconductor package **DECIPHER**, developed by Erik S. Wright and colleagues.

If you use FLIPPer with `--engine detectrepeats`, please cite:

Cho, S.-T. & Wright, E. S. (2025).  
*Accurate detection of tandem repeats exposes ubiquitous reuse of biological sequences.*  
Nucleic Acids Research, 53(17), gkaf866.  
https://doi.org/10.1093/nar/gkaf866

---

### Metapredict

FLIPPer also utilises **Metapredict**, an intrinsic disorder prediction tool developed by Ryan Emenecker, Daniel Griffith, and Alex Holehouse.

Please cite:

Emenecker, R. J., Griffith, D., & Holehouse, A. S. (2022).  
*Metapredict V2: An update to metapredict, a fast, accurate, and easy-to-use predictor of consensus disorder and structure.*  
bioRxiv.  
https://www.biorxiv.org/content/10.1101/2022.06.06.494887v1  

Emenecker, R. J., Griffith, D., & Holehouse, A. S. (2021).  
*Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure.*  
Biophysical Journal, 120(20), 4312–4319.  
https://www.cell.com/biophysj/fulltext/S0006-3495(21)00725-6  

---

## ⚠️ Important Notes

- The **FASTA record ID must not exceed 50 characters**.  
  Longer IDs will interfere with sequence filtering after repeat detection.
- **Migrating from the pre-`--engine` FLIPPer or from the old `FLIPPer-DetectRepeats/` fork**: both have been merged into this single `FLIPPer.py` with `--engine {xstream,detectrepeats}` (default `xstream`, so existing XSTREAM-based `--non-interactive` invocations keep working unchanged). One behavioural fix came along with the merge: the DetectRepeats fork's `--serine`/`--alanine` filters were previously near no-ops (a Biopython percent-scale bug meant almost every sequence passed regardless of threshold) - they now filter correctly under `--engine detectrepeats`, the same way they always have under `--engine xstream`, which may reduce candidate counts for existing DetectRepeats-based workflows.
- **`--non-interactive` is now a deprecated no-op**: running with no flags always behaved this way by default from v2.1 onwards. Existing scripts passing `--non-interactive` keep working unchanged; the flag just does nothing anymore. Scripts/HPC jobs that relied on the *old* bare `python3 FLIPPer.py` launching the interactive prompts should add `--interactive` to keep that behaviour.
