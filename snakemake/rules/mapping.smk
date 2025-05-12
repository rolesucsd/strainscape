"""
Mapping rules for the iHMP pipeline.
Handles read mapping, sorting, and filtering of BAM files.
"""

rule bwa_index:
    input: 
        fa = FILTERED_CONTIGS("{patient}")
    output: 
        amb = FILTERED_CONTIGS("{patient}") + ".amb",
        ann = FILTERED_CONTIGS("{patient}") + ".ann",
        bwt = FILTERED_CONTIGS("{patient}") + ".bwt",
        pac = FILTERED_CONTIGS("{patient}") + ".pac",
        sa = FILTERED_CONTIGS("{patient}") + ".sa"
    conda: CONDA_MAPPING
    shell: 
        "bwa index {input.fa}"

rule map_reads_bwa:
    input:
        ref = FILTERED_CONTIGS("{patient}"),
        amb = FILTERED_CONTIGS("{patient}") + ".amb",
        ann = FILTERED_CONTIGS("{patient}") + ".ann",
        bwt = FILTERED_CONTIGS("{patient}") + ".bwt",
        pac = FILTERED_CONTIGS("{patient}") + ".pac",
        sa = FILTERED_CONTIGS("{patient}") + ".sa",
        r1 = lambda wc: get_samples(wc.patient)[wc.sample]
    output:
        bam_all = MAP_BAM("{patient}", "{sample}"),
        sorted = SORT_BAM("{patient}", "{sample}"),
        flagstat = FLAGSTAT("{patient}", "{sample}"),
        index = SORT_BAM("{patient}", "{sample}") + ".bai"
    threads: 8
    conda: CONDA_MAPPING
    shell:
        r"""
        mkdir -p $(dirname {output.bam_all})
        bwa mem -t {threads} {input.ref} {input.r1} \
          | samtools view -b - > {output.bam_all}
        samtools flagstat {output.bam_all} > {output.flagstat}
        samtools view -b -F 4 -q 1 {output.bam_all} | samtools sort -o {output.sorted}
        samtools index {output.sorted}
        """ 