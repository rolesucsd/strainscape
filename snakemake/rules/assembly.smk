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

rule megahit_coassembly:
    """
    Perform co-assembly of multiple samples using MEGAHIT.
    
    This rule assembles reads from multiple samples into a single assembly
    using MEGAHIT. It uses the following parameters:
    - k-mer range: 27-77
    - k-mer step: 10
    - Merge level: 20,0.95
    - Minimum contig length: 200bp
    
    Input:
        reads: List of read files from all samples
    Output:
        fa: Final assembly file
    Resources:
        mem: 200GB
        cpu: 8 threads
    """
    input:
        lambda wc: list(get_samples(wc.patient).values())
    output: 
        fa = COASSEMBLY("{patient}")
    params:
        outdir = ASSEMBLY_DIR("{patient}"),
        kmin = 27,
        kmax = 77,
        kstep = 10,
        merge = "20,0.95",
        minlen = 200
    conda:
        config['conda_envs']['megahit']
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
        fna = BAKTA_FNA("{patient}")
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