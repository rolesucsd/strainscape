"""
MAGs rules for the iHMP pipeline.
Handles depth calculation and binning of metagenomic data.

This module contains rules for processing metagenomic assembled genomes (MAGs), including:
- Depth calculation using jgi_summarize_bam_contig_depths
- Binning using MetaBAT2
- Quality assessment of bins

Dependencies:
- Sorted BAM files from mapping
- Filtered contigs from assembly
- MetaBAT2 for binning
"""

import os

rule jgi_summarize_depths:
    """
    Calculate contig depths using jgi_summarize_bam_contig_depths.
    
    This rule calculates the depth of coverage for each contig across all
    samples. The depth information is used by MetaBAT2 for binning.
    
    Input:
        bams: List of sorted BAM files
    Output:
        depth: Depth file containing coverage information
    """
    input:
        bams = lambda wc: [SORT_BAM(wc.patient, s) for s in read_samples(wc.patient)]
    output:
        depth = DEPTH_FILE("{patient}")
    conda:
        config['conda_envs']['metabat2']
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams}
        """

rule metabat2_binning:
    """
    Perform binning using MetaBAT2.
    
    This rule uses MetaBAT2 to bin contigs into metagenomic assembled
    genomes (MAGs) based on:
    - Coverage depth
    - Tetranucleotide frequency
    - Contig length
    
    Parameters:
    - Minimum contig length: 1500bp
    - Maximum probability threshold: 95%
    - Minimum score: 60
    - Maximum edges: 200
    
    Input:
        contig: Filtered contigs file
        depth: Depth file from jgi_summarize_depths
    Output:
        directory: Directory containing binned contigs
    Resources:
        mem: 50GB
    """
    input:
        contig = FILTERED_CONTIGS("{patient}"),
        depth  = DEPTH_FILE("{patient}")
    output:
        directory(PATIENT_BIN_DIR("{patient}"))
    params:
        outdir = PATIENT_BIN_DIR("{patient}")
    conda:
        config['conda_envs']['metabat2']
    resources:
        mem = "50G"
    shell:
        r"""
        mkdir -p {params.outdir}
        metabat2 -i {input.contig} \
                 -a {input.depth} \
                 -o {params.outdir}/bin \
                 -s 1500 -m 1500 --maxP 95 --minS 60 --maxEdges 200 \
                 --seed 1 --saveCls
        """

rule checkm2:
    """
    Classify and assess MAGs using CheckM2.

    This rule runs CheckM2 to estimate the taxonomy and quality of bins.

    Input:
        bin_dir: Directory with binned MAGs (FASTA files)
    Output:
        results: Directory with CheckM2 output including taxonomy and QC
    """
    input:
        bin_dir = PATIENT_BIN_DIR("{patient}")
    output:
        results = directory(CHECKM2_OUT("{patient}")),
        tsv = os.path.join(CHECKM2_OUT("{patient}"), "quality_report.tsv")
    params:
        db=config['reference']['checkm2_db']
    conda:
        config['conda_envs']['checkm2']
    resources:
        mem = "50G"
    shell:
        """
        checkm2 predict -x fa --database_path {params.db} -i {input.bin_dir} -o {output.results}
        """

rule make_bin_txt:
    """
    Generate a bin.txt file (scaffold→bin mapping) and merge
    CheckM2 metrics (Completeness, Contamination, Genome_Size).
    """
    input:
        bin_dir    = PATIENT_BIN_DIR("{patient}"),
        checkm2_tsv = os.path.join(CHECKM2_OUT("{patient}"), "quality_report.tsv")
    output:
        bin_file = BIN_TXT("{patient}")
    log:
        os.path.join(config['paths']['log_dir'], "make_bin_txt_{patient}.log")
    shell:
        """
        python ../strainscape/make_bin_txt.py \
          --bin_dir {input.bin_dir} \
          --checkm2-tsv {input.checkm2_tsv} \
          --output_file {output.bin_file} \
          --log_file {log} 2>&1
        """
