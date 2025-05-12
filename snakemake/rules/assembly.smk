"""
Assembly rules for the iHMP pipeline.
Handles co-assembly, filtering, and annotation of metagenomic data.
"""

rule megahit_coassembly:
    input:
        lambda wc: list(get_samples(wc.patient).values())
    output: 
        fa = COASSEMBLY("{patient}")
    params:
        outdir = os.path.join(ASSEMBLY_DIR, "{patient}"),
        kmin = 27,
        kmax = 77,
        kstep = 10,
        merge = "20,0.95",
        minlen = 200
    conda: CONDA_MEGAHIT
    threads: 8
    resources: mem = "200G"
    shell:
        r"""
        mkdir -p {params.outdir}
        megahit \
          --force \
          -r $(printf "%s," {input}) \
          --min-contig-len {params.minlen} \
          --k-min {params.kmin} --k-max {params.kmax} --k-step {params.kstep} \
          --merge-level {params.merge} \
          -o {params.outdir}
        mv {params.outdir}/final.contigs.fa {output.fa}
        """

rule filter_contigs:
    input: 
        fa = COASSEMBLY("{patient}")
    output: 
        fa = FILTERED_CONTIGS("{patient}")
    conda: CONDA_SEQKIT
    shell: 
        "seqkit seq -m 1000 {input.fa} -o {output.fa}"

rule make_stb:
    input: 
        fa = FILTERED_CONTIGS("{patient}")
    output: 
        stb = STB_FILE("{patient}")
    run:
        with open(input.fa) as fin, open(output.stb, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    fout.write(f"{line[1:].split()[0]}\t{wildcards.patient}\n")

rule bakta_coassembly:
    input: 
        fa = FILTERED_CONTIGS("{patient}")
    output: 
        touch(os.path.join(REF_BAKTA_DIR, "{patient}", "bakta.done"))
    params:
        db = config["reference"]["bakta_db"],
        outdir = os.path.join(REF_BAKTA_DIR, "{patient}"),
        prefix = "{patient}"
    conda: CONDA_BAKTA
    resources: mem = "50G"
    shell:
        "mkdir -p {params.outdir} && bakta --db {params.db} --output {params.outdir} "
        "--prefix {params.prefix} --keep-contig-headers --force {input.fa}" 