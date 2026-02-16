# 🔬 Fast Linker Identification Pipeline for Pyrenoids (FLIPPer)

FLIPPer runs on **FASTA-formatted protein sequences** to identify candidate linker regions associated with pyrenoids.

All FASTA files must be placed in the **root directory** (the same directory as `FLIPPer.py`).  
When executed, FLIPPer scans the root directory for valid FASTA files and processes them using the defined pipeline settings.

---

## 🚀 Usage

```bash
python3 FLIPPer.py
```

Running this command from the terminal or command prompt launches the **parameterisation interface**.

---

## 📂 Outputs

For each analysed FASTA file, FLIPPer generates the following outputs:

### 1️⃣ Metapredict Plots
- PDF graphs containing:
  - Metapredict disorder scores  
  - Predicted AlphaFold2 pLDDT values  

---

### 2️⃣ XSTREAM Outputs
- Tandem repeat detection results
- Viewable as a linked HTML file:
  ```
  *_out_2.html
  ```

---

### 3️⃣ Candidate Sequences
- FASTA file containing all candidate linker sequences  
- CSV file containing the same sequences  

---

### 4️⃣ Variables File
- Text file containing the user-defined variables used during the pipeline run  

---

### 5️⃣ Sequence Analysis
- CSV files containing sequence analysis results generated during Step 1 of sequence filtering  

---

## 📦 Requirements

- Python ≥ 3.10  
- biopython  
- matplotlib  
- beautifulsoup4  
- pandas  
- metapredict  
- cython  

You can install dependencies using:

```bash
pip install biopython matplotlib beautifulsoup4 pandas metapredict cython
```

---

## 🧬 Third-Party Software & Citations

### XSTREAM

This pipeline utilises **XSTREAM**, a tandem repeat detection algorithm developed by Aaron Newman and James Cooper.

If you use FLIPPer, please cite:

Newman, A. M. & Cooper, J. B. (2007).  
*XSTREAM: a practical algorithm for identification and architecture modeling of tandem repeats in protein sequences.*  
BMC Bioinformatics, 8, 382.  
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-382

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
  Longer IDs will interfere with sequence filtering after XSTREAM repeat detection.