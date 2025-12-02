#!/usr/bin/awk -f

#-----------------------------------------------------------------------------------------------
#              FROM DELLY VCF, STRUCTURAL VARIANT (SV) DATA ANALYSIS
#
#
# OVERVIEW OF WHAT BEEN DONE:
# - PASS filtering applied
# - PRECISE filtering isn't applied only during child sample identification
# - Extracted CHROM, POS, END, SVTYPE, PE, SR, SV_LENGTH
# - Extracted genotypes from sample1, sample2, sample3
# - Child identification (s1, s2, s3) done using Mendelian violation.
# - De novo detection (child = sample3) done.
# - Parents were inferred (who father, who mother) using chrX heterozygosity
# - Child sex via chrX heterozygosity only
# - Variant statistics:
#    1. Count of variant types (SVTYPE)
#    2. Counts per chromosome
#    3. Genotype combination frequencies (sample1_sample2_sample3)
#    4. Bi-allelic vs multi-allelic ALT alleles
#    5. SNV vs SV counts
#-----------------------------------------------------------------------------------

BEGIN {
    FS  = "\t"
    OFS = "\t"

    print "CHROM","START","END","SV_TYPE","LENGTH_bp","Paired_end_PE","Split_end_SR",
          "sample1_HG00512","sample2_HG00513","sample3_HG00514" \
          > "output/SV_summary.txt"

    print "CHROM","START","END","SV_TYPE","LENGTH_bp","Paired_end_PE","Split_end_SR",
          "sample1_HG00512","sample2_HG00513","sample3_HG00514" \
          > "output/denovo_variants_precise.txt"
          
    print "CHROM","START","END","SV_TYPE","LENGTH_bp","PairedEnd_PE","SplitEnd_SR",
          "sample1_HG00512","sample2_HG00513","sample3_HG00514" \
          > "output/denovo_variants_imprecise.txt"
    
    print "=== Trio Summary===" \
          > "output/summary_stats.txt"

    print "" > "SV_summary.avinput"
    print "" > "denovo_variants_precise.avinput"
    print "" > "denovo_variants_imprecise.avinput"
    
    
    #counter for total variants
    total_variants = 0

    # child identification counters
    c1 = c2 = c3 = 0
    
    # chrX heterozygosity counters
    totalX = 0
    hetA = hetB = hetC = 0

    print "=== PROCESSING STARTED ==="
}


#==== SKIPPING THE VCF HEADERS===================
/^#/ { next }


#=========FILTERING PASS==========================
$7 != "PASS" { next }


#===========MAIN PROCESSING BLOCK============
{
    chrom = $1
    pos   = $2
    ref   = $4
    alter = $5
    info  = $8

#===========INITIALIZING VALUES==============

    svtype = "NA"
    endpos = pos
    pe = sr = 0
    precise = 0
    imprecise = 0
    
    
    
    # PARSING INFO
    n = split(info, element, ";")
    for (i=1; i<=n; i++) {
        if (element[i] == "PRECISE") precise = 1  
        if (element[i] == "IMPRECISE") imprecise = 1
        
        split(element[i], key_value, "=")
        key = key_value[1]; val = key_value[2]

        if (key=="SVTYPE") svtype = val
        else if (key=="END") endpos = val + 0
        else if (key=="PE") pe = val + 0
        else if (key=="SR") sr = val + 0
    }
    
    # DETERMINING LENGTH
    length_bp = endpos - pos
    
    #---------------------VARIANT READ COUNT SUMMARY---------------------------
    total_variants++
    if (precise == 1) precise_count++
    if (imprecise == 1) imprecise_count++

    #-----------------------GENOTYPES ---------------------------------------
    split($10, a, ":"); s1 = a[1]   # sample1
    split($11, b, ":"); s2 = b[1]   # sample2
    split($12, c, ":"); s3 = c[1]   # sample3

    if (s1=="./." || s2=="./." || s3=="./.") next


    #===========CHILD IDENTIFICATION (DEVIATIONS)==============
    
    # if Sample1 is child (deviations)
    if ((s1 == "1/1" && s2 == "0/0" && s3 == "0/0") || (s1 == "0/1" && s2 == "0/0" && s3 == "0/0"))c1++

    # if Sample2 is child (deviations)
    if ((s2 == "1/1" && s1 == "0/0" && s3 == "0/0") || (s2 == "0/1" && s1 == "0/0" && s3 == "0/0")) c2++

    # if Sample3 is child (deviations)
    if ((s3 == "1/1" && s1 == "0/0" && s2 == "0/0") || (s3 == "0/1" && s1 == "0/0" && s2 == "0/0")) c3++
    
  
    #===========LOADING INTO SUMMARY FILE==============
        print chrom, pos, endpos, svtype, length_bp, pe, sr, s1, s2, s3 \
              >> "output/SV_summary.txt"
          
        print chrom, pos, endpos, "N", "<"svtype">" \
              >> "SV_summary.avinput"

    
    #===========BI-ALLELIC AND MULTI-ALLELIC COUNT============
    #PRECISE FILTER APPLIED
   
    if (precise == 1) {
        # multi-allelic if ALT contains a comma
        if (index(alter, ",") > 0)
            multiallelic++
        else
            biallelic++
    
    }
    
    #===========SNV AND SV COUNT=================================
    #PRECISE FILTER APPLIED
    
    if (precise == 1) {
        # SNV count 
        if (alter !~ /,/ && length(ref)==1 && length(alter)==1)
            snv_count++

        # SV count 
        if (alter !~ /,/ && ((length(ref)>=1 && length(alter)>1) || (length(ref)>1 && length(alter)>=1) || (length(ref)>1 && length(alter)>1)))
        sv_count++

    }
        
    
    #===========DE NOVO VARIANTs=====================
    if ((s1=="0/0") && (s2=="0/0") && (s3=="0/1" || s3=="1/1")) {
    
   
    # PRECISE VARIANTS
    
        if (precise == 1) {
            print chrom, pos, endpos, svtype, length_bp, pe, sr, s1, s2, s3 \
                >> "output/denovo_variants_precise.txt"

            print chrom, pos, endpos, "N", "<"svtype">" \
                >> "denovo_variants_precise.avinput"
    }

   
    # IMPRECISE VARIANTS
    if (imprecise == 1 || precise == 1) {
        print chrom, pos, endpos, svtype, length_bp, pe, sr, s1, s2, s3 \
              >> "output/denovo_variants_imprecise.txt"

        print chrom, pos, endpos, "N", "<"svtype">" \
              >> "denovo_variants_imprecise.avinput"
    }
}
    
    #===========ChrX ANALYSIS FOR PARENT & SEX DETERMINATION===============
    # PRECISE FILTER APPLIED
    if (precise == 1) {
        if (chrom=="chrX") {
            totalX++
            if (s1=="0/1" || s1=="1/0") hetA++
            if (s2=="0/1" || s2=="1/0") hetB++
            if (s3=="0/1" || s3=="1/0") hetC++
        }
    }

    #===========STATS: VARIANT COUNTS PER TYPE, CHROMOSOME, GENOTYPE COMBINATION==============
    # PRECISE FILTER APPLIED
    if (precise == 1) {
        typeCount[svtype]++
        chromCount[chrom]++
        combination[s1 "_" s2 "_" s3]++
    }
}

#===========END BLOCK: PRINTING ALL SUMMARIES===================
END { 
        
    #===========CHILD IDENTIFICATION==============
    
    #for printing in the terminal
    print "Possible child = Sample1 (HG00512), deviations =", c1
    print "Possible child = Sample2 (HG00513), deviations =", c2
    print "Possible child = Sample3 (HG00514), deviations =", c3
    print "The sample with the least possible deviations is the child sample."

    #for the txt file
    print "sample1 deviations:", c1 >> "output/summary_stats.txt"
    print "sample2 deviations:", c2 >> "output/summary_stats.txt"
    print "sample3 deviations:", c3 >> "output/summary_stats.txt"
    print "The sample with the least possible deviations is the child sample." >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"

    
        
    #===========VARIANT READ COUNT SUMMARY======================
   
    print "\n=== Variant Read Count Summary ===" >> "output/summary_stats.txt"
    print "Total variant reads:", total_variants + 0 >> "output/summary_stats.txt"
    print "Precise reads:", precise_count + 0 >> "output/summary_stats.txt"
    print "Imprecise reads:", imprecise_count + 0 >> "output/summary_stats.txt"
    
    
    #===========NO. OF DIFFERENT TYPES OF ALLELES====================
    
    print "\n===Allele Structure Counts===" >> "output/summary_stats.txt"
    print "Bi-allelic variants:", biallelic + 0 >> "output/summary_stats.txt"
    print "Multi-allelic variants:", multiallelic + 0 >> "output/summary_stats.txt"


    #===========SNV AND SV====================    
    print "\n=== SNV / SV Counts===" >> "output/summary_stats.txt"
    print "SNVs (Single Nucleotide Variants):", snv_count + 0 >> "output/summary_stats.txt"
    print "SVs (Structural Variants):", sv_count + 0 >> "output/summary_stats.txt"
    
    
    #===========VARIANT TYPE COUNTS===============    
    print "\n===Variant Types Count===" >> "output/summary_stats.txt"
    for (t in typeCount)
        print t, typeCount[t] >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"

    
    #===========GENOTYPE COMBINATION COUNTS=====================
    print "\n===Genotype Combination Counts (Parent1_Parent2_Child)===" >> "output/summary_stats.txt"
    for (com in combination)
        print com, combination[com] >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"

    
    #===========VARIANTS PER CHROMOSOME========================
    print "\n===Variants per Chromosome===" >> "output/summary_stats.txt"
    for (chr in chromCount)
        print chr, chromCount[chr] >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"
    
    

    #===========ChrX SUMMARY=========================
    print "\n==chrX Heterozygosity Analysis==" >> "output/summary_stats.txt"
    print "Total chrX variants:", totalX >> "output/summary_stats.txt"
    print "ParentA (sample1, HG00512) het:", hetA >> "output/summary_stats.txt"
    print "ParentB (sample2, HG00513) het:", hetB >> "output/summary_stats.txt"
    print "Child   (sample3, HG00514) het:", hetC >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"
    
    #for printing in the terminal
    print "\n==chrX Heterozygosity Analysis==" 
    print "Total chrX variants:", totalX 
    print "ParentA (sample1, HG00512) het:", hetA
    print "ParentB (sample2, HG00513) het:", hetB
    print "Child   (sample3, HG00514) het:", hetC
 


    #===========INFERRING MOTHER / FATHER================
    if (hetA > hetB) {
        mother="sample1"; father="sample2"
    } else if (hetB > hetA) {
        mother="sample2"; father="sample1"
    } else {
        mother="Undetermined"; father="Undetermined"
    }

    print "\n===Inferred Parents===" >> "output/summary_stats.txt"
    print "Mother (likely):", mother >> "output/summary_stats.txt"
    print "Father (likely):", father >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"


    #for printing in the terminal
    print "\n===Inferred Parents===="
    print "Mother (likely):", mother 
    print "Father (likely):", father 
    
    
    
    
    #===========SEX OF THE CHILD===================================
    
    print "\n====Inferred Sex of Child (sample3)====" >> "output/summary_stats.txt"
    if (hetC > 0)
        print "Child = FEMALE (heterozygous chrX detected)" >> "output/summary_stats.txt"
    else
        print "Child = MALE (no chrX heterozygosity)" >> "output/summary_stats.txt"
    print "" >> "output/summary_stats.txt"

    
    #for printing in the terminal
    print "\n====Inferred Sex of Child (sample3)====" 
    if (hetC > 0)
        print "Child = FEMALE (heterozygous chrX detected)" 
    else
        print "Child = MALE (no chrX heterozygosity observed)" 
        
        
    #===========OUTPUT SUMMARY==================
    print "\n=== OUTPUT GENERATED ==="
    print "output/SV_summary.txt"
    print "output/denovo_variants_precise.txt"
    print "output/denovo_variants_imprecise.txt"
    print "SV_summary.avinput"
    print "denovo_variants_precise.avinput"
    print "denovo_variants_imprecise.avinput"
    print "output/summary_stats.txt"
}
