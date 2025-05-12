"""
Quality control rules for the iHMP pipeline.
Handles read quality assessment and filtering.
"""

rule run_fastqc:
    input:
        r1 = lambda wc: get_samples(wc.patient)[wc.sample]
    output:
        html = FASTQC_HTML("{patient}", "{sample}"),
        zip = FASTQC_ZIP("{patient}", "{sample}")
    threads: 8
    conda: CONDA_QC
    shell:
        r"""
        fastqc {input.r1} -o $(dirname {output.html}) -t {threads}
        """

rule run_multiqc:
    input:
        fastqc = expand(FASTQC_ZIP("{patient}", "{sample}"), 
                       sample=lambda wc: get_samples(wc.patient).keys())
    output:
        report = MULTIQC_REPORT("{patient}")
    threads: 8
    conda: CONDA_QC
    shell:
        r"""
        multiqc {input.fastqc} -o {output.report}
        """

rule run_trim_galore:
    input:
        r1 = lambda wc: get_samples(wc.patient)[wc.sample]
    output:
        r1_trimmed = TRIMMED_READS("{patient}", "{sample}")
    threads: 8
    conda: CONDA_QC
    shell:
        r"""
        trim_galore --quality 20 --length 50 --output_dir $(dirname {output.r1_trimmed}) \
          --cores {threads} {input.r1}
        """ 