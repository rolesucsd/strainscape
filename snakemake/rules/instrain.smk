"""
InStrain rules for the iHMP pipeline.
Handles strain-level analysis using InStrain.
"""

rule run_instrain_profile:
    input:
        fa = FILTERED_CONTIGS("{patient}"),
        bam = SORT_BAM("{patient}", "{sample}")
    output:
        isdb = INSTRAIN_DB("{patient}", "{sample}")
    threads: 8
    conda: CONDA_INSTRAIN
    shell:
        r"""
        inStrain profile {input.bam} {input.fa} -o {output.isdb} -p {threads} \
          --database_mode
        """

rule run_instrain_compare:
    input:
        isdb = expand(INSTRAIN_DB("{patient}", "{sample}"), 
                     sample=lambda wc: get_samples(wc.patient).keys())
    output:
        comparison = INSTRAIN_COMPARE("{patient}")
    threads: 8
    conda: CONDA_INSTRAIN
    shell:
        r"""
        inStrain compare -i {input.isdb} -o {output.comparison} -p {threads}
        """

rule run_instrain_snvs:
    input:
        isdb = INSTRAIN_DB("{patient}", "{sample}")
    output:
        snvs = INSTRAIN_SNVS("{patient}", "{sample}")
    threads: 8
    conda: CONDA_INSTRAIN
    shell:
        r"""
        inStrain snv_call -i {input.isdb} -o {output.snvs} -p {threads}
        """ 