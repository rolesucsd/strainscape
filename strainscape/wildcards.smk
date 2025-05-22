"""
Wildcard function definitions for the iHMP pipeline.
These functions define the file paths used throughout the workflow.
DATA_DIR is set from the config file, so all outputs can be redirected by the user.
"""

import os
import yaml

def get_data_dir():
    """
    Get the data directory from the config file.
    Returns the data_dir path from the config file.
    """
    return config["paths"]["data_dir"]

# Assembly wildcards
def COASSEMBLY(patient):
    return f"{get_data_dir()}/assembly/{patient}/final_contigs.fa"

def FILTERED_CONTIGS(patient):
    return f"{get_data_dir()}/assembly/{patient}/contigs_gt1kb.fa"

def STB_FILE(patient):
    return f"{get_data_dir()}/assembly/{patient}/assembly.stb"

def ASSEMBLY_DIR(patient):
    return f"{get_data_dir()}/assembly/{patient}"

# Mapping wildcards
def MAP_BAM(patient, sample):
    return f"{get_data_dir()}/mapping/{patient}/{sample}.all.bam"

def SORT_BAM(patient, sample):
    return f"{get_data_dir()}/mapping/{patient}/{sample}.filtered.sorted.bam"

def FLAGSTAT(patient, sample):
    return f"{get_data_dir()}/mapping/{patient}/{sample}.filtered.flagstat.txt"

def COMBINED_FLAGSTAT():
    return f"{get_data_dir()}/mapping/combined/combined.flagstat.txt"

# Depth and binning wildcards
def DEPTH_FILE(patient):
    return f"{get_data_dir()}/depth/{patient}.depth.txt"

def PATIENT_BIN_DIR(patient):
    return f"{get_data_dir()}/bins/{patient}"

def BIN_TXT(patient):
    return f"{get_data_dir()}/assembly/{patient}/bins.txt"

def CHECKM2_OUT(patient):
    return f"{get_data_dir()}/checkm2/{patient}"

# InStrain wildcards
def INSTR_DIR(patient):
    return f"{get_data_dir()}/instrain/{patient}"

def SNV_FILE(patient, sample):
    return f"{get_data_dir()}/instrain/{patient}/each/{sample}/output/{sample}_SNVs.tsv"

def SCF_FILE(patient, sample):
    return f"{get_data_dir()}/instrain/{patient}/each/{sample}/output/{sample}_scaffold_info.tsv"

def COMBINED_SCAFFOLD_INFO(patient):
    return f"{get_data_dir()}/instrain/{patient}/combined_scaffold_info.tsv"

def COMBINED_SNV_INFO(patient):
    return f"{get_data_dir()}/instrain/{patient}/combined_SNV_info.tsv"

# Annotation wildcards
def REF_BAKTA_DIR(patient):
    return f"{get_data_dir()}/bakta/{patient}"

def BAKTA_TSV(patient):
    return f"{get_data_dir()}/bakta/{patient}/{patient}.tsv"

def BAKTA_FNA(patient):
    return f"{get_data_dir()}/bakta/{patient}/{patient}.fna"

def MAGS_BIN(patient, bin_id):
    return f"{get_data_dir()}/bins/{patient}/{bin_id}.fa"

# SNP tracking wildcards
def PROCESSED_SCAFFOLD(patient):
    return f"{get_data_dir()}/instrain/{patient}/processed_scaffolds.txt"

def FILTERED_MUTATIONS(patient):
    return f"{get_data_dir()}/instrain/{patient}/SNV_filtered.txt"

def MUTATION_TRENDS(patient):
    return f"{get_data_dir()}/instrain/{patient}/SNV_filtered_trend.txt"

def MUTATION_TRAJECTORIES(patient):
    return f"{get_data_dir()}/instrain/{patient}/SNV_filtered_trend_trajectory.txt"

def MAPPED_MUTATIONS(patient):
    return f"{get_data_dir()}/instrain/{patient}/SNV_filtered_trend_trajectory_mapped.txt"

def MUTATION_TYPES(patient):
    return f"{get_data_dir()}/instrain/{patient}/SNV_filtered_trend_trajectory_mapped_types.txt"

def STYLE(patient):
    return f"{get_data_dir()}/assembly/{patient}/style_never_made.joke"