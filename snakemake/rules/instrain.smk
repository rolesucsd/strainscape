"""

This module contains rules for performing strain-level analysis using InStrain,
including:
- Running InStrain profile on individual samples
- Combining scaffold information across samples
- Combining SNV information across samples
- Processing scaffold information for downstream analysis

Dependencies:
- Sorted BAM files from mapping
- Filtered contigs from assembly
- STB files for strain tracking
"""

# Import wildcard functions
from scripts.wildcards import (
    FILTERED_CONTIGS, SORT_BAM, STB_FILE,
    SNV_FILE, SCF_FILE, INSTR_DIR,
    COMBINED_SCAFFOLD_INFO, COMBINED_SNV_INFO,
    PROCESSED_SCAFFOLD
)

rule inStrain_profile:
    """
    Run InStrain profile on individual samples.
    
    This rule performs strain-level analysis on each sample using InStrain,
    generating SNV and scaffold information files. It uses the following
    parameters:
    - Minimum coverage: 10x
    - Database mode enabled
    - All reads considered for pairing
    
    Input:
        bam_sorted: Sorted BAM file from mapping
    Output:
        snv: SNV information file
        scaffold: Scaffold information file
    Resources:
        mem: 100GB
        cpu: 8 threads
    """
    input:
        bam_sorted = SORT_BAM("{patient}", "{sample}")
    output:
        snv = SNV_FILE("{patient}", "{sample}"),
        scaffold = SCF_FILE("{patient}", "{sample}")
    params:
        outdir = os.path.join(INSTR_DIR("{patient}"), "each", "{sample}", "output"),
        genome = FILTERED_CONTIGS("{patient}"),
        stb = STB_FILE("{patient}")
    conda:
        config['conda_envs']['instrain']
    threads: 8
    resources: mem = "100G"
    shell: 
        r"""
        mkdir -p {params.outdir}
        inStrain profile {input.bam_sorted} {params.genome} \
            -s {params.stb} -o {params.outdir} \
            --min_cov 10 -p {threads} --database_mode -d \
            --pairing_filter all_reads \
            --skip_plot_generation
        """

rule combine_scaffold_info:
    """
    Combine scaffold information across samples.
    
    This rule merges scaffold information from all samples for a patient,
    creating a combined file for downstream analysis. 

    Input:
        scaffold_files: Individual scaffold information files
    Output:
        combined: Combined scaffold information file
    """
    input:
        lambda wc: [SCF_FILE(wc.patient, s) for s in get_samples(wc.patient)]
    output:
        combined = COMBINED_SCAFFOLD_INFO("{patient}")
    run:
        import pandas as pd, os
        dfs = [
            pd.read_csv(f, sep="\t").assign(Sample=os.path.basename(os.path.dirname(f)))
            for f in input if os.path.getsize(f)
        ]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined, sep="\t", index=False)

rule combine_SNV_info:
    """
    Combine SNV information across samples.
    
    This rule merges SNV information from all samples for a patient,
    creating a combined file for downstream analysis. 
    
    Input:
        snv_files: Individual SNV information files
    Output:
        combined: Combined SNV information file
    """
    input:
        lambda wc: [SNV_FILE(wc.patient, s) for s in get_samples(wc.patient)]
    output:
        combined = COMBINED_SNV_INFO("{patient}")
    run:
        import pandas as pd, os
        dfs = [
            pd.read_csv(f, sep="\t").assign(Sample=os.path.basename(os.path.dirname(f)))
            for f in input if os.path.getsize(f)
        ]
        pd.concat(dfs, ignore_index=True).to_csv(output.combined, sep="\t", index=False)

rule process_scaffolds:
    """
    Process combined scaffold information for downstream analysis.
    
    This rule processes the combined scaffold information by:
    1. Filtering scaffolds based on length, coverage, and breadth thresholds
    2. Merging with STB and metadata information
    3. Calculating summary statistics
    
    Input:
        scaffold_file: Combined scaffold information file
        stb_file: STB file for bin mapping
        metadata_file: Metadata file with sample information
    
    Output:
        processed_file: Processed scaffold information file
    
    Parameters:
        min_length: Minimum scaffold length (default: 1000)
        min_coverage: Minimum coverage threshold (default: 10.0)
        min_breadth: Minimum breadth threshold (default: 0.4)
        threads: Number of threads to use (default: 4)
        chunksize: Size of chunks for processing (default: 10000)
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        scaffold_file = rules.combine_scaffold_info.output.combined_file,
        stb_file = rules.make_stb.output.stb_file,
        metadata_file = config['metadata']['file']
    output:
        processed_file = os.path.join(config['paths']['results'], "processed_scaffolds.tsv")
    params:
        min_length = config.get("min_scaffold_length", 1000),
        min_coverage = config.get("min_scaffold_coverage", 5.0),
        min_breadth = config.get("min_scaffold_breadth", 0.4),
        threads = config.get("threads", 4),
        chunksize = config.get("chunksize", 10000)
    resources:
        mem_mb = 16000,
        threads = 4
    log:
        os.path.join(config['paths']['log_dir'], "process_scaffolds.log")
    conda:
        config['conda_envs']['instrain']
    shell:
        """
        python {workflow.basedir}/scripts/process_scaffolds.py \
            --scaffold_file {input.scaffold_file} \
            --stb_file {input.stb_file} \
            --metadata_file {input.metadata_file} \
            --output_file {output.processed_file} \
            --min_length {params.min_length} \
            --min_coverage {params.min_coverage} \
            --min_breadth {params.min_breadth} \
            --threads {params.threads} \
            --chunksize {params.chunksize} \
            --log_file {log} 2>&1
        """

rule merge_snv_info:
    """
    Merge SNV information with STB, metadata, and filtered scaffolds.
    
    This rule combines SNV information with bin mapping (STB), patient metadata,
    and filtered scaffold information to ensure SNVs only map to scaffolds that
    passed previous filtering criteria.
    
    Input:
        snv_file: Combined SNV information file
        stb_file: STB file for bin mapping
        metadata_file: Metadata file with sample information
        processed_scaffold_file: Processed scaffold file containing filtered scaffolds
    
    Output:
        merged_file: Merged SNV information file
    
    Parameters:
        threads: Number of threads to use (default: 4)
        chunksize: Size of chunks for processing (default: 10000)
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        snv_file = rules.combine_SNV_info.output.combined_file,
        stb_file = rules.make_stb.output.stb_file,
        metadata_file = config['metadata']['file'],
        processed_scaffold_file = rules.process_scaffolds.output.processed_file
    output:
        merged_file = os.path.join(config['paths']['results'], "merged_snvs.tsv")
    params:
        threads = config.get("threads", 4),
        chunksize = config.get("chunksize", 10000)
    resources:
        mem_mb = 16000,
        threads = 4
    log:
        os.path.join(config['paths']['log_dir'], "merge_snv_info.log")
    conda:
        config['conda_envs']['instrain']
    shell:
        """
        python {workflow.basedir}/scripts/filter_mutations.py \
            --snv_file {input.snv_file} \
            --stb_file {input.stb_file} \
            --metadata_file {input.metadata_file} \
            --processed_scaffold_file {input.processed_scaffold_file} \
            --output_file {output.merged_file} \
            --threads {params.threads} \
            --chunksize {params.chunksize} \
            --log_file {log} 2>&1
        """
