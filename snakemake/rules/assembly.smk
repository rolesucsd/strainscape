"""
Assembly rules for the iHMP pipeline.
Handles co-assembly, filtering, and annotation of metagenomic data.

This module contains rules for processing metagenomic assemblies, including:
- Co-assembly of multiple samples using MEGAHIT
- Contig filtering based on length
- STB file generation for strain tracking
- Gene annotation using Bakta

Dependencies:
- Raw sequencing reads
- MEGAHIT for assembly
- SeqKit for contig filtering
- Bakta for gene annotation
"""

rule concat_reads_for_coassembly:
    """
    Concatenate all R1s and all R2s for a patient into two gz files.
    Inputs are assumed to be *.fastq.gz; gzip multi-member streams are OK.
    """
    input:
        r1 = lambda wc: [read_samples(wc.patient)[s]["R1"] for s in sorted(read_samples(wc.patient))],
        r2 = lambda wc: [read_samples(wc.patient)[s]["R2"] for s in sorted(read_samples(wc.patient))]
    output:
        r1 = COASSEMBLY_R1("{patient}"),
        r2 = COASSEMBLY_R2("{patient}")
    threads: 1
    shell:
        r"""
        mkdir -p "$(dirname {output.r1})"
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """

rule metaspades_coassembly:
    """
    Co-assemble all timepoints for this patient into a single assembly.

    Choices we made:
    - --meta: metagenomic mode.
    - -k 33..127: broad k-mer sweep improves contiguity for mixed-abundance communities.
    - --only-assembler: skip read error correction to reduce compute; assumes upstream QC.
      (If you don't do quality trimming/adapter removal earlier, drop this flag.)
    - --phred-offset 33: modern Illumina default.
    - --cov-cutoff auto: SPAdes heuristic to drop ultra-low cov junk.
    - --min-contig: enforce MIN_CONTIG_LEN here to reduce downstream binning noise.
    """
    input:
        r1 = COASSEMBLY_R1("{patient}"),
        r2 = COASSEMBLY_R2("{patient}")
    output:
        COASSEMBLY("{patient}")
    params:
        outdir = COASSEMBLY_DIR("{patient}"),
        kmer   = "33,55,77,99"
    conda:
        config['conda_envs']['metaspades']
    threads: 16
    shell:
        r"""
        metaspades.py \
          -1 {input.r1} \
          -2 {input.r2} \
          --meta \
          --only-assembler \
          -k {params.kmer} \
          --phred-offset 33 \
          -o "{params.outdir}"
        """

rule filter_contigs:
    """
    Filter contigs based on length.
    
    This rule removes short contigs from the assembly, keeping only those
    longer than 1000bp. This helps reduce noise and improve downstream
    analysis quality.
    
    Input:
        fa: Raw assembly file
    Output:
        fa: Filtered assembly file
    """
    input: 
        fa = COASSEMBLY("{patient}")
    output: 
        fa = FILTERED_CONTIGS("{patient}")
    conda:
        config['conda_envs']['seqkit']
    shell: 
        "seqkit seq -m 1000 {input.fa} -o {output.fa}"

rule make_stb:
    """
    Generate STB file for strain tracking.
    
    This rule creates a scaffold-to-bin (STB) file that maps contigs to
    their source patient. This information is used by InStrain for strain
    tracking and analysis.
    
    Input:
        fa: Filtered assembly file
    Output:
        stb: STB file
    """
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
    """
    Annotate assembly using Bakta.
    
    This rule performs gene annotation on the filtered assembly using Bakta.
    It generates both a TSV file with gene annotations and a FNA file with
    the annotated sequences.
    
    Input:
        fa: Filtered assembly file
    Output:
        tsv: Bakta annotation file
        fna: Bakta annotated sequences file
    Resources:
        mem: 50GB
    """
    input: 
        fa = FILTERED_CONTIGS("{patient}")
    output: 
        tsv = BAKTA_TSV("{patient}"),
        fna = BAKTA_FNA("{patient}"),
        faa = BAKTA_FAA("{patient}")
    params:
        db = config["reference"]["bakta_db"],
        outdir = REF_BAKTA_DIR("{patient}"),
        prefix = "{patient}"
    conda:
        config['conda_envs']['bakta']
    resources: mem = "50G"
    shell:
        """
        mkdir -p {params.outdir} && bakta --db {params.db} --output {params.outdir} \
            --prefix {params.prefix} --keep-contig-headers --force {input.fa}
        """

rule bakta_prefix_faa:
    input:
        faa = BAKTA_FAA("{patient}")
    output:
        pref = BAKTA_FAA_PREFIX("{patient}")
    shell:
        r"""
        awk 'BEGIN{{OFS=""}}
             /^>/{{sub(/^>/,">" "{wildcards.patient}" "|"); print; next}}
             {{print}}' {input.faa} > {output.pref}
        """

rule concat_bakta_proteins_all:
    input:
        expand(BAKTA_FAA_PREFIX("{patient}"),
               patient=PATIENT_IDS)
    output:
        faa = os.path.join(COMBINED_MMSEQ2_DIR(), "combined_proteins.faa")
    shell:
        "mkdir -p results/mmseqs && cat {input} > {output.faa}"

rule mmseqs_cluster:
    input:
        faa = os.path.join(COMBINED_MMSEQ2_DIR(), "combined_proteins.faa")
    output:
        db  = os.path.join(COMBINED_MMSEQ2_DIR(), "db"),
        clu = os.path.join(COMBINED_MMSEQ2_DIR(), "clu"),
        tsv = os.path.join(COMBINED_MMSEQ2_DIR(), "clusters.tsv"),
        rep = os.path.join(COMBINED_MMSEQ2_DIR(), "cluster_reps.faa")
    conda: 
        config['conda_envs']['mmseq2']
    threads: 16
    shell:
        r"""
        mmseqs createdb {input.faa} {output.db}
        mmseqs cluster {output.db} {output.clu} results/mmseqs/tmp \
            --min-seq-id 0.6 --cov-mode 1 -c 0.8 --threads {threads}
        mmseqs createtsv {output.db} {output.db} {output.clu} {output.tsv}
        mmseqs createseqfiledb {output.db} {output.clu} results/mmseqs/rep_db
        mmseqs result2flat {output.db} results/mmseqs/rep_db results/mmseqs/rep_db {output.rep}
        """

rule eggnog_map_cluster_reps:
    input:
        rep = os.path.join(COMBINED_MMSEQ2_DIR(), "cluster_reps.faa")
    output:
        tsv = os.path.join(COMBINED_MMSEQ2_DIR(), "eggnog_annotation.tsv"),
    params:
        db = config["reference"]["eggnog_db"]
    conda:
        config['conda_envs']['eggnog']
    threads: 16
    shell:
        r"""
        mkdir -p results/eggnog
        emapper.py -i {input.rep} -o results/eggnog/eggnog \
          --output_dir results/eggnog --data_dir {params.db} \
          --cpu {threads} --override --decorate_gff no --itype proteins
        mv results/eggnog/eggnog.emapper.annotations {output.tsv}
        """
