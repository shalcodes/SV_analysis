import pandas as pd
import subprocess
import os


def annotate_sv(
    input_file,
    output_file,
    gene_bed="hg38_refGene.bed",
    exon_bed="hg38_exons.bed",
    clinvar_bed="clinvar_SV.bed"):
    
# for generating clinvar_SV.bed, we used https://www.ncbi.nlm.nih.gov/clinvar/?term=%22structural+variant%22 
# (downloaded the txt file from here)

    """
    This annotation function is using the following:
       - Gene overlaps
       - Exon overlaps
       - Checking priority (Exonic > Intronic > Intergenic)
       - ClinVar pathogenicity category

    Expected input_file format (ANNOVAR input format, AVINPUT):
        chrom   start   end   ref   alt

    ClinVar BED format:
        chrom   start   end   germline_classification
        
    ClinVar BED format (with condition):
        chrom   start   end   germline_classification  condition
    """

# ----------------------------------------------------------------
# 1. If file is empty, then will produce an empty annotated file
# ----------------------------------------------------------------
    if os.path.getsize(input_file) == 0:
        print(f"{input_file} is empty — writing empty annotation file.")

        empty_cols = [
            "chrom","start","end","ref","alt",
            "Function","Gene","Priority",
            "clinvar_germline_classification",
            "clinvar_condition"
            
        ]
        pd.DataFrame(columns=empty_cols).to_csv(output_file, index=False)
        return pd.DataFrame()


# ------------------------------------------------
# 2. Loading input SVs (chrom, start, end, ref, alt)
# ------------------------------------------------
    df = pd.read_csv(
        input_file, sep="\t", header=None,
        names=["chrom","start","end","ref","alt"]
    )

    # Replaces missing REF or ALT
    df["ref"] = df["ref"].fillna("N")
    df["alt"] = df["alt"].fillna("<NA>")
    

    # Creates BED file for BEDTools (internally only)
    df["start0"] = df["start"] - 1  
    # convert for BEDTools. internal only, not in output.
    
    sv_bed = input_file + ".tmp.bed"
    df[["chrom","start0","end","alt"]].to_csv(
        sv_bed, sep="\t", header=False, index=False
    )  
    
    
    
# ----------------------------------------------------
# 3. Gene overlaps
# ----------------------------------------------------
    gene_overlap_file = input_file + ".gene_overlap"
    subprocess.run(
        ["bedtools", "intersect", "-a", sv_bed, "-b", gene_bed, "-wa", "-wb"],
        stdout=open(gene_overlap_file, "w"),
        check=False
    )

    gene_cols = [
        "sv_chrom","sv_start0","sv_end","sv_alt",
        "gene_chrom","gene_start","gene_end","gene_name"
    ]

    try:
        gene_df = pd.read_csv(gene_overlap_file, sep="\t", names=gene_cols)
    except pd.errors.EmptyDataError:
        gene_df = pd.DataFrame(columns=gene_cols)

# ----------------------------------------------------
# 4. Exon overlaps
# ----------------------------------------------------
    exon_overlap_file = input_file + ".exon_overlap"
    subprocess.run(
        ["bedtools", "intersect", "-a", sv_bed, "-b", exon_bed, "-wa", "-wb"],
        stdout=open(exon_overlap_file, "w"),
        check=False
    )

    exon_cols = [
        "sv_chrom","sv_start0","sv_end","sv_alt",
        "exon_chrom","exon_start","exon_end","exon_gene"
    ]

    try:
        exon_df = pd.read_csv(exon_overlap_file, sep="\t", names=exon_cols)
    except pd.errors.EmptyDataError:
        exon_df = pd.DataFrame(columns=exon_cols)

    exonic_sv = set(zip(
        exon_df["sv_chrom"],
        exon_df["sv_start0"],
        exon_df["sv_end"]
    ))
    
    
    
# ----------------------------------------------------
# 5. Functional classification
# ----------------------------------------------------
    def classify(row):
        key = (row.chrom, row.start0, row.end)

        if key in exonic_sv:
            return "exonic"

        hits = gene_df[
            (gene_df["sv_chrom"] == row.chrom) &
            (gene_df["sv_start0"] == row.start0) &
            (gene_df["sv_end"] == row.end)
        ]

        if not hits.empty:
            return "intronic"

        return "intergenic"

    df["Function"] = df.apply(classify, axis=1)
    
    

# ----------------------------------------------------
# 6. Assigning overlapping gene(s)
# ----------------------------------------------------
    def get_gene(row):
        hits = gene_df[
            (gene_df["sv_chrom"] == row.chrom) &
            (gene_df["sv_start0"] == row.start0) &
            (gene_df["sv_end"] == row.end)
        ]

        if hits.empty:
            return "-"

        unique_genes = sorted(hits["gene_name"].unique())
        return ";".join(unique_genes)

    df["Gene"] = df.apply(get_gene, axis=1)
        
    

# ----------------------------------------------------
# 7. Priority score
# ----------------------------------------------------
    df["Priority"] = df["Function"].map({
        "exonic": "High",
        "intronic": "Medium",
        "intergenic": "Low"
    })


# ----------------------------------------------------
# 8. ClinVar pathogenicity annotation
# ----------------------------------------------------
    clin_file = input_file + ".clinvar_overlap"

    subprocess.run([
        "bedtools", "intersect",
        "-a", sv_bed,
        "-b", clinvar_bed,
        "-wa", "-wb"
    ], stdout=open(clin_file, "w"), check=False)

    clin_cols = [
        "sv_chrom","sv_start0","sv_end","sv_alt",
        "c_chrom","c_start","c_end","c_germ"
    ]

    try:
        clin_df = pd.read_csv(clin_file, sep="\t", names=clin_cols)
    except pd.errors.EmptyDataError:
        clin_df = pd.DataFrame(columns=clin_cols)

    def lookup_clinvar(row):
        hits = clin_df[
            (clin_df["sv_chrom"] == row.chrom) &
            (clin_df["sv_start0"] == row.start0) &
            (clin_df["sv_end"] == row.end)
        ]

        if hits.empty:
            return "-"

        return ";".join(sorted(hits["c_germ"].unique()))

    df["clinvar_germline_classification"] = df.apply(lookup_clinvar, axis=1)
    
    
# ----------------------------------------------------
# 8. ClinVar Clinical Condition annotation
# ----------------------------------------------------
    clin_cond_file = input_file + ".clinvar_condition_overlap"

    subprocess.run([
        "bedtools", "intersect",
        "-a", sv_bed,
        "-b", "clinvar_SV_condition.bed",
        "-wa", "-wb"
    ], stdout=open(clin_cond_file, "w"), check=False)

    cond_cols = [
        "sv_chrom","sv_start0","sv_end","sv_alt",
        "c_chrom","c_start","c_end","condition","GermlineClass"
    ]

    try:
        cond_df = pd.read_csv(clin_cond_file, sep="\t", names=cond_cols)
    except pd.errors.EmptyDataError:
        cond_df = pd.DataFrame(columns=cond_cols)

    def lookup_condition(row):
        hits = cond_df[
            (cond_df["sv_chrom"] == row.chrom) &
            (cond_df["sv_start0"] == row.start0) &
            (cond_df["sv_end"] == row.end)
        ]
        if hits.empty:
            return "-"
        return ";".join(sorted(hits["condition"].unique()))

    df["clinvar_condition"] = df.apply(lookup_condition, axis=1)

    
    
# ----------------------------------------------------
# 9. Final output
# ----------------------------------------------------
    final_cols = [
        "chrom","start","end","ref","alt",
        "Function","Gene","Priority",
        "clinvar_germline_classification", 
        "clinvar_condition"
    ]

    df[final_cols].to_csv(output_file, index=False)

    print(f"Annotation written to {output_file}")
    return df[final_cols]