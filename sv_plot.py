#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from matplotlib import ticker

# -------------------------------------------------------------------
# PREPARE OUTPUT FOLDER
# -------------------------------------------------------------------
PLOT_DIR = "plots"
os.makedirs(PLOT_DIR, exist_ok=True)

sns.set(style="whitegrid", font_scale=1.2)


# -------------------------------------------------------------------
# 1. SAVING FIGURES
# -------------------------------------------------------------------
def save_plot(fig, name):
    path = os.path.join(PLOT_DIR, name)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"{path} is saved.")
    plt.close(fig)


def add_value_labels(ax, fmt="%.0f"):
    """
    Add value labels on top of each bar in a barplot / countplot.
    Ensures counts appear as whole numbers.
    """
    for container in ax.containers:
        ax.bar_label(container, fmt=fmt)


# -------------------------------------------------------------------
# 2. PLOTING FOR SV ANNOTATED FILES
# -------------------------------------------------------------------
def plot_sv_annotation(csv_file, tag="SV_summary"):
    print(f"==={csv_file} is being loaded===")
    df = pd.read_csv(csv_file)

    # Basic columns
    df["SV_size"] = df["end"] - df["start"]

    # -----------------------------------------------------------------
    # PLOT 1: Variant Function Distribution
    # -----------------------------------------------------------------
    if tag != "SV_exonic":
        fig, ax = plt.subplots(figsize=(8, 5))
        order = df["Function"].value_counts().index
        sns.countplot(
            data=df,
            x="Function",
            order=order,
            ax=ax,
            palette="Blues"
        )
        ax.set_title(f"{tag}: Variant Classification by Genomic Function")
        ax.set_ylabel("Number of Variants")
        ax.set_xlabel("Function Category")
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        add_value_labels(ax)
        save_plot(fig, f"{tag}_function_distribution.png")
    else:
        print(f"---Skipping Function Distribution plot for {tag} (all exonic)---")
    

    # ---------------------------------------------------------------------
    # PLOT 2: Chromosome-wise SV Distribution
    # ---------------------------------------------------------------------
    def chrom_sort_key(ch):
        name = ch.replace("chr", "")
        if name.isdigit():
            return (0, int(name))
        if name == "X":
            return (1, 23)
        if name == "Y":
            return (1, 24)
        return (2, name)

    chrom_order = sorted(df["chrom"].unique(), key=chrom_sort_key)

    fig, ax = plt.subplots(figsize=(9, 5))
    sns.countplot(
        data=df,
        x="chrom",
        order=chrom_order,
        ax=ax,
        palette="GnBu"
    )
    ax.set_title(f"{tag}: Variants per Chromosome")
    ax.set_ylabel("Count")
    ax.set_xlabel("Chromosome")
    plt.xticks(rotation=45)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    add_value_labels(ax)
    save_plot(fig, f"{tag}_chromosome_distribution.png")

    # -----------------------------------------------------------------
    # PLOT 3: Top Most Affected Genes (genes on X, counts on Y)
    # -----------------------------------------------------------------
    df_genes = df[df["Gene"] != "-"]

    if not df_genes.empty:
        exploded = df_genes["Gene"].str.split(";", expand=True).stack()
        top_genes = exploded.value_counts().nlargest(20)

        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(
            x=top_genes.index,
            y=top_genes.values,
            ax=ax,
            palette="Oranges_r"
        )
        ax.set_title(f"{tag}: Top Affected Genes")
        ax.set_xlabel("Gene")
        ax.set_ylabel("Variant Count")
        plt.xticks(rotation=60, ha="right")
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        add_value_labels(ax)
        save_plot(fig, f"{tag}_top_genes.png")

    print(f"-----SV Plots for {csv_file} is done----")
    

    
# -------------------------------------------------------------------
# 3. SOME PLOTS FROM summary_stats.txt 
# -------------------------------------------------------------------

def plot_summary_stats(stats_file="summary_stats.txt"):
    print(f"----Reading summary stats: {stats_file}----")
    with open(stats_file) as f:
        content = f.read()
        
    # Extract bi-allelic and multi-allelic
    bi = re.search(r"Bi-allelic variants:\s+(\d+)", content)
    multi = re.search(r"Multi-allelic variants:\s+(\d+)", content)
    
    # Extract structural variants and single nucleotide variants
    snv = re.search(r"SNVs.*:\s+(\d+)", content)
    sv = re.search(r"SVs.*:\s+(\d+)", content)
    

   
    #------------------------SNV vs SV-----------------------------
    
    if snv and sv:
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.barplot(
            x=["SNVs", "SVs"],
            y=[int(snv.group(1)), int(sv.group(1))],
            palette="Purples",
            ax=ax
        )
        ax.set_title("SNV vs SV Count")
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        add_value_labels(ax)
        save_plot(fig, "summary_snvs_vs_svs.png")
        
        
    
    # ----------------------Bi-allelic & MULTI-ALLELIC----------------------
    
    if bi and multi:
        fig, ax = plt.subplots(figsize=(7, 5))
        y_vals = [int(bi.group(1)), int(multi.group(1))]
        sns.barplot(
            x=["Bi-allelic", "Multi-allelic"],
            y=y_vals,
            ax=ax,
            palette="Reds"
        )
        ax.set_title("Allele Structure Count")
        ax.set_ylabel("Count")
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        add_value_labels(ax)
        save_plot(fig, "summary_allele_structure.png")
        
        
#-------------------PRECISE and IMPRECISE READS---------------------------

def plot_precise_imprecise_pie(stats_file="summary_stats.txt"):
    print(f"----Reading precise/imprecise counts from: {stats_file}----")

    with open(stats_file) as f:
        content = f.read()

    # Extract counts using regex
    precise = re.search(r"Precise reads:\s+(\d+)", content)
    imprecise = re.search(r"Imprecise reads:\s+(\d+)", content)

    if not precise or not imprecise:
        print("----Could not find precise/imprecise counts in summary_stats.txt----")
        return

    precise_count = int(precise.group(1))
    imprecise_count = int(imprecise.group(1))

    labels = ["Precise Reads", "Imprecise Reads"]
    values = [precise_count, imprecise_count]

    fig, ax = plt.subplots(figsize=(7, 7))

    wedges, texts, autotexts = ax.pie(
        values,
        labels=None,
        autopct="%1.1f%%",
        startangle=140,
        textprops={"color": "black"},
    )

    plt.setp(autotexts, size=14, weight="bold")
    ax.set_title("Precise vs Imprecise Variant Reads", fontsize=16)
    
    # ------- Adding legend with counts ---------
    ax.legend(
        wedges,
        labels,
        title="Variant Read Type",
        loc="center left",
        bbox_to_anchor=(1.05, 0.5),
        fontsize=12
    )

    save_plot(fig, "summary_precise_imprecise_pie.png")

    
# -------------------------------------------------------------------
# 4. ALT Variant Type Count (DEL / INS / DUP) plots
# -------------------------------------------------------------------
def plot_alt_counts(csv_file, tag="ALT_counts"):
    print(f"=== Loading ALT counts from {csv_file} ===")
    df = pd.read_csv(csv_file)

    alt_counts = df["alt"].value_counts()
    
    rename_map = {
        "<DEL>": "DEL",
        "<INS>": "INS",
        "<DUP>": "DUP"
    }
    labels = [rename_map.get(x, x) for x in alt_counts.index]
    values = alt_counts.values

    # -------------------- PLOT --------------------
    fig, ax = plt.subplots(figsize=(8, 5))
    sns.barplot(
        x=labels,
        y=values,
        ax=ax,
        palette="Oranges_r"   
    )
    
    if tag == "SV_exonic":
        ax.set_title("Counts of SV types in Annotated Exonic Variants")
        
    elif tag == "SV_pathlink":
        ax.set_title("Counts of SV types in Annotated Variants with Pathogenic Links")
        
    elif tag == "SV_annotated":
        ax.set_title("Counts of SV types in All Annotated Variants")
        
    ax.set_xlabel("Variant Types")
    ax.set_ylabel("Count")

    plt.xticks(rotation=0)

    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    add_value_labels(ax)

    plt.tight_layout()
    save_plot(fig, f"{tag}_alt_type_counts.png")
    print("=== ALT count plot saved ===")
  
    
# -------------------------------------------------------------------
# 4. MAIN FUNCTION 
# -------------------------------------------------------------------
def run_all_plots():

    print("\n=== Plotting Full SV Summary ===")
    if os.path.exists("SV_summary_annotated.csv"):
        plot_sv_annotation("SV_summary_annotated.csv", tag="SV_summary")

    print("\n=== Plotting Imprecise de novo ===")
    if os.path.exists("denovo_variants_imprecise_annotated.csv"):
        plot_sv_annotation("denovo_variants_imprecise_annotated.csv", tag="denovo_imprecise")

    print("\n=== Plotting Some Info from Summary Stats ===")
    if os.path.exists("summary_stats.txt"):
        plot_summary_stats("summary_stats.txt")

    print("\n=== Plotting EXONIC variants ===")
    if os.path.exists("SV_summary_annotated_exonic.csv"):
        plot_sv_annotation("SV_summary_annotated_exonic.csv", tag="SV_exonic")

    print("\n=== Plotting PATHOGENIC variants ===")
    if os.path.exists("SV_summary_annotated_pathLink.csv"):
        plot_sv_annotation("SV_summary_annotated_pathLink.csv", tag="SV_pathogenic")
        
    print("\n=== Plotting Precise vs Imprecise Pie Chart ===")
    if os.path.exists("summary_stats.txt"):
        plot_precise_imprecise_pie("summary_stats.txt")
    
    print("\n=== Plotting ALT Value Counts ===")
    if os.path.exists("SV_summary_annotated_exonic.csv"):
        plot_alt_counts("SV_summary_annotated_exonic.csv", tag="SV_exonic")
        
    if os.path.exists("SV_summary_annotated_pathLink.csv"):
        plot_alt_counts("SV_summary_annotated_pathLink.csv", tag="SV_pathlink")
        
    if os.path.exists("SV_summary_annotated.csv"):
        plot_alt_counts("SV_summary_annotated.csv", tag="SV_annotated")
    

if __name__ == "__main__":
    run_all_plots()
    
    


