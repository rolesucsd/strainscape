"""
Mapping rules for the iHMP pipeline.
Handles read mapping, sorting, and filtering of BAM files.

This module contains rules for processing sequencing reads, including:
- BWA indexing of reference genomes
- Read mapping using BWA-MEM
- BAM file sorting and filtering
- Flagstat generation for mapping statistics
- Combining flagstat files across samples

Dependencies:
- Filtered contigs from assembly
- Raw sequencing reads
- BWA and SAMtools tools
"""

rule bwa_index:
    """
    Create BWA index for reference genome.
    
    This rule generates the necessary index files for BWA-MEM alignment.
    
    Input:
        fa: Filtered contigs file
    Output:
        amb: Ambiguity index file
        ann: Annotation index file
        bwt: BWT index file
        pac: Packed sequence index file
        sa: Suffix array index file
    """
    input: 
        fa = FILTERED_CONTIGS("{patient}")
    output: 
        amb = FILTERED_CONTIGS("{patient}") + ".amb",
        ann = FILTERED_CONTIGS("{patient}") + ".ann",
        bwt = FILTERED_CONTIGS("{patient}") + ".bwt",
        pac = FILTERED_CONTIGS("{patient}") + ".pac",
        sa = FILTERED_CONTIGS("{patient}") + ".sa"
    conda:
        config['conda_envs']['mapping']
    shell: 
        "bwa index {input.fa}"

rule map_reads_bwa:
    """
    Map reads to reference genome using BWA-MEM.
    
    This rule performs the following steps:
    1. Maps reads using BWA-MEM
    2. Converts SAM to BAM
    3. Filters unmapped reads and low-quality mappings
    4. Sorts BAM file
    5. Generates flagstat statistics
    
    Input:
        ref: Reference genome
        amb, ann, bwt, pac, sa: BWA index files
        r1: Forward reads
    Output:
        bam_all: Unfiltered BAM file
        sorted: Sorted and filtered BAM file
        flagstat: Mapping statistics
        index: BAM index file
    Resources:
        cpu: 8 threads
    """
    input:
        ref = FILTERED_CONTIGS("{patient}"),
        amb = FILTERED_CONTIGS("{patient}") + ".amb",
        ann = FILTERED_CONTIGS("{patient}") + ".ann",
        bwt = FILTERED_CONTIGS("{patient}") + ".bwt",
        pac = FILTERED_CONTIGS("{patient}") + ".pac",
        sa = FILTERED_CONTIGS("{patient}") + ".sa",
        r1 = lambda wc: read_samples(wc.patient)[wc.sample]
    output:
        bam_all = MAP_BAM("{patient}", "{sample}"),
        sorted = SORT_BAM("{patient}", "{sample}"),
        flagstat = FLAGSTAT("{patient}", "{sample}"),
        index = SORT_BAM("{patient}", "{sample}") + ".bai"
    threads: 8
    conda:
        config['conda_envs']['mapping']
    shell:
        r"""
        mkdir -p $(dirname {output.bam_all})
        bwa mem -t {threads} {input.ref} {input.r1} \
          | samtools view -b - > {output.bam_all}
        samtools flagstat {output.bam_all} > {output.flagstat}
        samtools view -b -F 4 -q 1 {output.bam_all} | samtools sort -o {output.sorted}
        samtools index {output.sorted}
        """

rule combine_flagstats:
    """
    Combine flagstat files from mapping outputs.
    
    This rule aggregates mapping statistics from all samples into a single
    file for easier analysis and reporting. It processes the following
    statistics:
    - Total reads
    - Mapped reads
    - Mapping quality
    - Duplicate reads
    
    Input:
        flagstats: Individual flagstat files
    Output:
        combined: Combined flagstat file
    """
    input:
        flagstats = expand(FLAGSTAT("{patient}", "{sample}"),
                         patient=config["patient_ids"],
                         sample=lambda wc: read_samples(wc.patient).keys())
    output:
        combined = COMBINED_FLAGSTAT()
    params:
        base_dir = "results/mapping"
    conda:
        config['conda_envs']['mapping']
    script:
        "../scripts/combine_flagstats.py" 