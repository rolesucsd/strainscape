"""
Wildcard function definitions for the iHMP pipeline.
These functions define the file paths used throughout the workflow.
"""

# Assembly wildcards
def COASSEMBLY(patient):
    return f"data/assembly/{patient}/final_contigs.fa"

def FILTERED_CONTIGS(patient):
    return f"data/assembly/{patient}/contigs_gt1kb.fa"

def STB_FILE(patient):
    return f"data/assembly/{patient}/assembly.stb"

def ASSEMBLY_DIR(patient):
    return f"data/assembly/{patient}"

# Mapping wildcards
def MAP_BAM(patient, sample):
    return f"data/mapping/{patient}/{sample}.all.bam"

def SORT_BAM(patient, sample):
    return f"data/mapping/{patient}/{sample}.filtered.sorted.bam"

def FLAGSTAT(patient, sample):
    return f"data/mapping/{patient}/{sample}.filtered.flagstat.txt"

def COMBINED_FLAGSTAT():
    return "data/mapping/combined/combined.flagstat.txt"

# Depth and binning wildcards
def DEPTH_FILE(patient):
    return f"data/depth/{patient}.depth.txt"

def PATIENT_BIN_DIR(patient):
    return f"data/bins/{patient}"

# InStrain wildcards
def INSTR_DIR(patient):
    return f"data/instrain/{patient}"

def SNV_FILE(patient, sample):
    return f"data/instrain/{patient}/each/{sample}/output/{sample}_SNVs.tsv"

def SCF_FILE(patient, sample):
    return f"data/instrain/{patient}/each/{sample}/output/{sample}_scaffold_info.tsv"

def COMBINED_SCAFFOLD_INFO(patient):
    return f"data/instrain/{patient}/combined/scaffold_info.tsv"

def COMBINED_SNV_INFO(patient):
    return f"data/instrain/{patient}/combined/snv_info.tsv"

# Annotation wildcards
def REF_BAKTA_DIR(patient):
    return f"data/bakta/{patient}"

def BAKTA_TSV(patient):
    return f"data/bakta/{patient}/{patient}.tsv"

def MAGS_BIN(patient, bin_id):
    return f"data/bins/{patient}/{bin_id}.fa"

# SNP tracking wildcards
def FILTERED_MUTATIONS(patient):
    return f"data/instrain/{patient}/SNV_filtered.txt"

def MUTATION_TRENDS(patient):
    return f"data/instrain/{patient}/SNV_filtered_trend.txt"

def MUTATION_TRAJECTORIES(patient):
    return f"data/instrain/{patient}/SNV_filtered_trend_trajectory.txt"

def MAPPED_MUTATIONS(patient):
    return f"data/instrain/{patient}/SNV_filtered_trend_trajectory_cluster_mapped.txt"

def MUTATION_TYPES(patient):
    return f"data/instrain/{patient}/SNV_filtered_trend_trajectory_cluster_mapped_mutation.txt" 