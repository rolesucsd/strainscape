"""
Annotation rules for the iHMP pipeline.
Handles functional and taxonomic annotation of genomes.
"""

rule run_prokka:
    input:
        fa = MAGS_DEREP("{patient}")
    output:
        gff = PROKKA_GFF("{patient}"),
        faa = PROKKA_FAA("{patient}"),
        ffn = PROKKA_FFN("{patient}")
    threads: 8
    conda: CONDA_ANNOTATION
    shell:
        r"""
        prokka --outdir $(dirname {output.gff}) --prefix $(basename {output.gff} .gff) \
          --cpus {threads} {input.fa}
        """

rule run_eggnog:
    input:
        faa = PROKKA_FAA("{patient}")
    output:
        emapper = EGGNOG("{patient}")
    threads: 8
    conda: CONDA_ANNOTATION
    shell:
        r"""
        emapper.py -i {input.faa} --output_dir {output.emapper} --cpu {threads}
        """

rule run_gt_db_tk:
    input:
        fa = MAGS_DEREP("{patient}")
    output:
        taxonomy = GTDB_TAXONOMY("{patient}")
    threads: 8
    conda: CONDA_ANNOTATION
    shell:
        r"""
        gtdbtk classify_wf --genome_dir $(dirname {input.fa}) \
          --out_dir {output.taxonomy} --cpus {threads}
        """ 