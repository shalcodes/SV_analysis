# Structural Variant (SV) Analysis Pipeline  
### *MSc Bioinformatics Project*

This repository contains a complete, end‑to‑end structural variant (SV) analysis workflow developed as part of an MSc Bioinformatics project.  
The pipeline processes a **Delly-generated VCF**, performs **sample detection (child and parents)**, **parent detection (who father and who mother)**, **sex determination of the child**,  **SV extraction**, **de novo detection**, **family‑based inference**, **BED-based annotation**, and **high‑quality visualisation**.

It combines shell scripting, AWK, Python, ANNOVAR database processing, and BEDTools-based comparisons.

---

# 🔬 Project Overview

This project demonstrates how to:

- Parse and process **structural variants** from a raw **Delly VCF**
- Extract PASS‑filtered variants and compute SV statistics
- Identify the **child sample** in a trio using Mendelian logic
- Detect **de novo structural variants**
- Infer **parents** and **child sex** using chrX heterozygosity
- Convert ANNOVAR reference files into BED files for annotation
- Annotate SVs using **gene**, **exon**, and **ClinVar SV** datasets
- Generate publication‑quality **plots** and **summary reports**

---

# 📁 Repository Structure

```
├── pipeline.py                 # Main script (orchestrates full workflow)
├── scripts/
│   ├── script.sh               # AWK-based DELLY VCF parser (SV extraction)
│   ├── annotate_sv.py          # Annotation using BEDTools & Python
│   └── bed.sh                  # Convert ANNOVAR refGene + ClinVar SV → BED
│
├── data/
│   ├── annovar/humandb/        # ANN OVAR reference files (not included here)
│   ├── bed_files/              # Generated BEDs (hg38_refGene.bed, exons, ClinVar)
│   └── DellyVariation.vcf      # Input VCF (example)
│
├── output/
│   ├── processed/              # Generated .txt and .avinput files
│   ├── stats/                  # Summary statistics
│   └── figures/                # Plots
│
├── requirements.txt
└── README.md
```

Users may delete all contents inside **output/** and regenerate them when running the workflow.

---

# 🧬 Methodology (High‑level)

The workflow contains **four main components**:

## **1. DELLY VCF Parsing (AWK — script.sh)**  
The VCF file is processed using an AWK script:

- Filters out non‑PASS variants  
- Extracts: `CHROM`, `POS`, `END`, `SVTYPE`, `PE`, `SR`, `SV length`  
- Extracts genotypes for three samples (HG00512, HG00513, HG00514)  
- Detects the **child** using Mendelian rules  
- Detects **de novo SVs**  
- Computes:  
  - Variant counts  
  - SV type distribution  
  - Chromosome enrichment  
  - SNV vs SV ratio  
  - Bi‑allelic vs multi‑allelic  
- Infers parents using **chrX heterozygosity**  
- Infers sex of child (presence/absence of chrX heterozygosity)  
- Generates:  
```
output/SV_summary.txt  
output/denovo_variants_precise.txt  
output/denovo_variants_imprecise.txt  
*.avinput  
output/summary_stats.txt
```

---

## **2. Reference BED Generation (bed.sh)**  
To enable annotation, the script converts reference datasets into BED format.

### **Datasets used**

| Source | File | Purpose |
|-------|------|---------|
| **ANNOVAR human genome database** | `annovar/humandb/hg38_refGene.txt` | Generate gene BED + exon BED |
| **NCBI ClinVar structural variant dataset** | `ClinVar_SV.txt` | Generate ClinVar SV BEDs |

These files are not uploaded here due to size limits, but the scripts expect them to be present in:

```
data/annovar/humandb/
```

### **Generated BED files**

| BED File | Description |
|---------|-------------|
| `hg38_refGene.bed` | Gene-level regions |
| `hg38_exons.bed` | Expanded exon-level regions |
| `clinvar_SV.bed` | ClinVar SVs filtered by variant type |
| `clinvar_SV_condition.bed` | ClinVar SVs annotated with condition + germline class |

These are later used with BEDTools to annotate DELLY variants.

---

## **3. Annotation (Python — annotate_sv.py)**  
The annotation script performs:

- **BEDTools intersect** between AVINPUT files and BED files
- Mapping SVs to:
  - Gene names  
  - Exon locations  
  - Known pathogenic/benign ClinVar regions  
- Generation of annotation tables (`*_annotated.csv`)

---

## **4. Visualisation (Python — sv_plot.py)**  
The script generates high‑quality figures:

- SV type barplots  
- Chromosomal density maps  
- PE/SR support analysis  
- De novo vs inherited patterns  
- Combination genotype matrix  
- Summary figures for reports/publications  

Saved under `output/figures/`.

---

# 📌 Input Files Used in This Project

### **1. DellyVariation.vcf**
A DELLY‑generated trio VCF used for demonstrating the entire workflow.

### **2. ANNOVAR Reference Database Files**
Downloaded from:  
https://annovar.openbioinformatics.org/en/latest/

Specifically used:  
- `hg38_refGene.txt`

### **3. ClinVar Structural Variant Dataset**
Downloaded manually from NCBI:  
https://www.ncbi.nlm.nih.gov/clinvar/?term=%22structural+variant%22

Saved as:  
```
annovar/humandb/ClinVar_SV.txt
```

---

# 🚀 Running the Full Pipeline

## **1. Install Dependencies**
```
pip install -r requirements.txt
```

Optional (recommended):
```
conda env create -f environment.yml
conda activate sv-pipeline
```

---

## **2. Make Scripts Executable**
```
chmod +x scripts/script.sh
chmod +x scripts/bed.sh
```

---

## **3. Generate BED Reference Files**
```
bash scripts/bed.sh
```

---

## **4. Run the Main Pipeline**
The entire workflow is controlled through:

```
python pipeline.py
```

This will:

- Parse VCF  
- Generate summary + .avinput files  
- Annotate variants  
- Create visuals  
- Save results into the **output/** directory  

---

# 🗂 Outputs Produced

| Folder | Contents |
|--------|----------|
| `output/processed/` | Parsed SV tables, AVINPUT files |
| `output/stats/` | Summary reports |
| `output/figures/` | Plots |
| root output | Denovo, Summary TXT files |

Users may delete all contents inside this folder and regenerate them.

---

# 🧠 Skills Demonstrated

- Next‑generation sequencing (NGS) analysis  
- Structural variant interpretation  
- AWK scripting  
- Python workflow automation  
- BEDTools‑based annotation  
- Data cleaning of genomic reference databases  
- Trio‑aware variant analysis  
- Visualisation and report generation  

---

# 📜 License

This repository is released under the **MIT License**.

---

# 👤 Contact  
**Shalini Majumder**
*sxm2220@student.bham.ac.uk / shalinimajumder24@gmail.com*
*MSc Bioinformatics, University of Birmingham*

