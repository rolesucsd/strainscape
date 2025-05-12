"""
MAGs (Metagenome-Assembled Genomes) rules for the iHMP pipeline.
Handles binning, refinement, and quality assessment of MAGs.
"""

rule run_metabat2:
    input:
        fa = FILTERED_CONTIGS("{patient}"),
        bam = expand(SORT_BAM("{patient}", "{sample}"), sample=lambda wc: get_samples(wc.patient).keys())
    output:
        bins = MAGS("{patient}")
    threads: 8
    conda: CONDA_MAGS
    shell:
        r"""
        metabat2 -i {input.fa} -a {input.bam} -o {output.bins} -t {threads}
        """

rule run_checkm:
    input:
        bins = MAGS("{patient}")
    output:
        quality = MAGS_QUALITY("{patient}")
    threads: 8
    conda: CONDA_MAGS
    shell:
        r"""
        checkm lineage_wf -t {threads} -x fa {input.bins} {output.quality}
        """

rule run_drep:
    input:
        bins = MAGS("{patient}")
    output:
        dereplicated = MAGS_DEREP("{patient}")
    threads: 8
    conda: CONDA_MAGS
    shell:
        r"""
        dRep dereplicate {output.dereplicated} -g {input.bins}/*.fa -p {threads}
        """ 