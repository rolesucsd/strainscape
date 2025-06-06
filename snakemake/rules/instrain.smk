"""

This module contains rules for performing strain-level analysis using InStrain,
including:
- Running InStrain profile on individual samples
- Combining scaffold information across samples
- Combining SNV information across samples
- Processing scaffold information for downstream analysis

Dependencies:
- Sorted BAM files from mapping
- Filtered contigs from assembly
- STB files for strain tracking
"""

rule inStrain_profile:
    """
    Run InStrain profile on individual samples.
    
    This rule performs strain-level analysis on each sample using InStrain,
    generating SNV and scaffold information files. It uses the following
    parameters:
    - Minimum coverage: 10x
    - Database mode enabled
    - All reads considered for pairing
    
    Input:
        bam_sorted: Sorted BAM file from mapping
    Output:
        snv: SNV information file
        scaffold: Scaffold information file
    Resources:
        mem: 100GB
        cpu: 8 threads
    """
    input:
        bam_sorted = SORT_BAM("{patient}", "{sample}")
    output:
        snv = SNV_FILE("{patient}", "{sample}"),
        scaffold = SCF_FILE("{patient}", "{sample}")
    params:
        outdir = os.path.join(INSTR_DIR("{patient}"), "each", "{sample}", "output"),
        genome = FILTERED_CONTIGS("{patient}"),
        stb = STB_FILE("{patient}")
    conda:
        config['conda_envs']['instrain']
    threads: 8
    resources: mem = "100G"
    shell: 
        r"""
        mkdir -p {params.outdir}
        inStrain profile {input.bam_sorted} {params.genome} \
            -s {params.stb} -o {params.outdir} \
            --min_cov 10 -p {threads} --database_mode -d \
            --pairing_filter all_reads \
            --skip_plot_generation
        """

rule combine_scaffold_info:
    """
    Combine scaffold information across samples.
    
    This rule merges scaffold information from all samples for a patient,
    creating a combined file for downstream analysis. 

    Input:
        scaffold_files: Individual scaffold information files
    Output:
        combined: Combined scaffold information file
    """
    input:
        lambda wc: [SCF_FILE(wc.patient, s) for s in read_samples(wc.patient)]
    output:
        combined = COMBINED_SCAFFOLD_INFO("{patient}")
    run:
        import os, sys
        import pandas as pd

        dfs = []
        expected_cols = None

        for f in input:
            # 1) skip zero-byte files
            if not os.path.getsize(f):
                print(f"[combine_scaffold_info] WARNING: {f} is zero-byte – skipped",
                      file=sys.stderr)
                continue

            # 2) try reading header first time to capture expected columns
            if expected_cols is None:
                try:
                    expected_cols = pd.read_csv(f, sep="\t", nrows=0).columns.tolist()
                except pd.errors.EmptyDataError:
                    # header-only or malformed—will be caught below
                    expected_cols = []

            # 3) read full table, catch truly empty‐content files
            try:
                df = pd.read_csv(f, sep="\t")
            except pd.errors.EmptyDataError:
                print(f"[combine_scaffold_info] WARNING: {f} has no columns – skipped",
                      file=sys.stderr)
                continue

            # 4) record sample name and keep
            df["Sample"] = os.path.basename(os.path.dirname(os.path.dirname(f)))
            dfs.append(df)

        # 5) write concatenated or fallback empty
        if dfs:
            pd.concat(dfs, ignore_index=True) \
              .to_csv(output.combined, sep="\t", index=False)
        else:
            # build empty DataFrame with the right schema
            cols = (expected_cols or []) + ["Sample"]
            pd.DataFrame(columns=cols) \
              .to_csv(output.combined, sep="\t", index=False)
            print(f"[combine_scaffold_info] WARNING: all input files empty → "
                  f"wrote empty {output.combined}",
                  file=sys.stderr)

rule combine_SNV_info:
    """
    Combine SNV information across samples, but don’t blow up
    if some of the per‐sample files are empty or unparseable.
    """
    input:
        lambda wc: [SNV_FILE(wc.patient, s) for s in read_samples(wc.patient)]
    output:
        combined = COMBINED_SNV_INFO("{patient}")
    run:
        import os, sys, pandas as pd

        dfs = []
        for fn in input:
            # skip zero‐byte files immediately
            if os.path.getsize(fn) == 0:
                print(f"[combine_SNV_info] WARNING: {fn} is zero‐byte, skipping", file=sys.stderr)
                continue
            # try to read; catch truly empty/no-column cases
            try:
                df = pd.read_csv(fn, sep="\t")
            except pd.errors.EmptyDataError:
                print(f"[combine_SNV_info] WARNING: {fn} has no columns or is unparseable, skipping", file=sys.stderr)
                continue

            # tag and collect
            df["Sample"] = os.path.basename(os.path.dirname(os.path.dirname(fn)))
            dfs.append(df)

        # write out
        if dfs:
            pd.concat(dfs, ignore_index=True) \
              .to_csv(output.combined, sep="\t", index=False)
        else:
            # emit a header‐only file so downstream rules still see a TSV
            pd.DataFrame(columns=["Sample"]) \
              .to_csv(output.combined, sep="\t", index=False)
            print(f"[combine_SNV_info] WARNING: all inputs empty → wrote header‐only {output.combined}",
                  file=sys.stderr)
