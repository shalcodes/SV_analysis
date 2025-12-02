# Structural Variant (SV) Analysis Pipeline
### *VCF → AWK Parsing → De Novo Detection → BED/ClinVar Annotation → Python Visualisation*

This repository contains a complete, reproducible workflow for **structural variant (SV) analysis** from a **Delly-generated VCF**, including:  
✔ SV extraction using AWK  
✔ De novo variant detection  
✔ ChrX-based parent inference & sex determination  
✔ ANNOVAR-style `.avinput` generation  
✔ Gene/exon/ClinVar annotation  
✔ Summary statistics & high-quality plots

---

## 📁 Repository Structure

```
├── script.sh              # AWK-based parser for DELLY VCF
├── annotate_sv.py         # annotates SVs using BED files
├── annovar
│   ├── humandb/             
│          ├── clinvar_SV.txt
│          ├── hg38_refGene.txt
│          └── hg38_refGeneMrna.fa
├── pipeline.py            # main pipeline controller
├── sv_plot.py             # custom plotting utilities
│
├── bed_files/             # gene/exon/clinvar annotations
├── example_vcf/           # optional
│
├── output/
│   ├── processed/             # parsed txt + avinput files
│   ├── stats/                 # summary statistics
│   └── figures/               # generated plots
│
├── requirements.txt
└── README.md
```

---

## 🧬 Overview of the Workflow

### **1. AWK Parsing of DELLY VCF**
The pipeline begins with an AWK script that:

- Applies **PASS** filtering  
- Extracts **CHROM, POS, END, SVTYPE, REF/ALT, PE, SR**  
- Computes **SV length**  
- Extracts **genotypes** for a trio (HG00512, HG00513, HG00514)  
- Performs **Mendelian violation checks** to identify the child  
- Detects **de novo variants**  
- Infers **mother/father** from chrX heterozygosity  
- Infers **child sex**  
- Generates:

```
output/SV_summary.txt
output/denovo_variants_precise.txt
output/denovo_variants_imprecise.txt
SV_summary.avinput
denovo_variants_precise.avinput
denovo_variants_imprecise.avinput
output/summary_stats.txt
```

---

## **2. Annotation (Python + BEDTools)**

`annotate_sv.py` annotates `.avinput` files with:
- Gene overlaps  
- Exonic vs intronic regions  
- ClinVar pathogenicity labels  
- Functional consequences  

Outputs:  
```
*_annotated.csv
```

---

## **3. Visualisation (`sv_plot.py`)**

The pipeline generates plots for:

- SV type distribution  
- Chromosomal density  
- Read support (PE/SR)  
- De novo vs inherited SVs  
- Allele structure (bi-allelic vs multi-allelic)  
- SNV vs SV counts  

Saved under:
```
output/figures/
```

---

## 🚀 Running the Pipeline

### **1. Install Dependencies**

```
pip install -r requirements.txt
```

Or with conda:

```
conda env create -f environment.yml
conda activate sv-pipeline
```

---

### **2. Ensure AWK is executable**

```
chmod +x scripts/script.sh
```

---

### **3. Run the Pipeline**

```
python src/pipeline.py
```

---

## 🔧 Requirements

- Python 3.8+
- AWK (GNU awk recommended)
- BEDTools
- ANNOVAR (optional)

Python dependencies (in `requirements.txt`):  
```
pandas
matplotlib
seaborn
```

---

## 📊 Example Results

- Counts of SV types  
- Genotype combination matrix  
- De novo variant lists  
- Chromosomal distribution plots  
- Parent & sex inference  
- ClinVar-annotated SV tables  

---

## 📜 Citation / Attribution

This project includes an AWK-based VCF parsing logic developed by **Shalini Majumder**.

---

## 📬 Contact

**Shalini Majumder**  
*MSc Bioinformatics, University of Birmingham*
