import subprocess
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import sv_plot as svp


#===========CONFIGURATION=============

# Delly VCF parsing script (.sh wrapper that calls awk)
DELLY_SCRIPT = "script.sh"            

# Input VCF
VCF_FILE = "DellyVariation.vcf"

# AVINPUT files produced by AWK
ALL_AVINPUT = "SV_summary.avinput"
DENOVO_AVINPUT_PRECISE = "denovo_variants_precise.avinput"
DENOVO_AVINPUT_IMPRECISE = "denovo_variants_imprecise.avinput"


#=============Python to Terminal====================

def run_cmd(cmd, label):
    print(f"\n=== {label} ===")
    try:
        subprocess.run(cmd, check=True)
        print(f"{label} completed.")
    except subprocess.CalledProcessError:
        print(f"{label} FAILED.")
        sys.exit(1)


#==============1. RUNNING AWK DELLY-PIPELINE=========================

run_cmd(
    ["awk", "-f", DELLY_SCRIPT, VCF_FILE],
    "===Running AWK==="
)

#=====================2. ANNOTATING USING BEDTOOLS==========================

from annotate_sv import annotate_sv

print("\n===Annotating SV_summary.avinput using BEDTools===")
annotate_sv("SV_summary.avinput", "output/SV_summary_annotated.csv")

print("\n===Annotating denovo_variants_precise.avinput using BEDTools===")
annotate_sv("denovo_variants_precise.avinput", "output/denovo_variants_precise_annotated.csv")

print("\n===Annotating denovo_variants_imprecise.avinput using BEDTools===")
annotate_sv("denovo_variants_imprecise.avinput", "output/denovo_variants_imprecise_annotated.csv")



#=================3. EXTRACTING EXONIC and PATHOGENIC SVs=====================

print("\n=== Extracting EXONIC variants ===")
sv = pd.read_csv("output/SV_summary_annotated.csv")
sv_exonic = sv[sv["Function"] == "exonic"]
sv_exonic.to_csv("output/SV_summary_annotated_exonic.csv", index=False)
print("output/SV_summary_annotated_exonic.csv is saved.")

print("\n=== Extracting Pathogenic / Likely Pathogenic variants ===")
sv_path = sv[
    sv["clinvar_germline_classification"].str.contains(
        "pathogenic|likely pathogenic",
        case=False,
        na=False
    )
]
sv_path.to_csv("output/SV_summary_annotated_pathLink.csv", index=False)
print("output/SV_summary_annotated_pathLink.csv is saved.")


#=======================4. RUNNING PLOTS==========================

print("\n=== Generating Plots ===")
svp.run_all_plots()

print("\n=== PIPELINE COMPLETED WITH PLOTS ===")



#===========DONE==================

print("\n============================================================")
print("FULL PIPELINE COMPLETED")
print("============================================================")
print("Generated files:")
print("1. SV_summary.txt")
print("2. SV_summary.avinput")
print("3. denovo_variants_precise.txt")
print("4. denovo_variants_imprecise.txt")
print("5. denovo_variants_precise.avinput")
print("6. denovo_variants_imprecise.avinput")
print("7. output/SV_summary_annotated.csv")
print("8. output/denovo_variants_precise_annotated.csv")
print("9. output/denovo_variants_imprecise_annotated.csv")
print("10. output/SV_summary_annotated_exonic.csv")
print("11. output/SV_summary_annotated_pathLink.csv")
print("============================================================")
