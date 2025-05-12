"""
Rules for SNP tracking and analysis in the iHMP pipeline.
"""

rule filter_mutations:
    """
    Filter mutations based on coverage and frequency thresholds.
    """
    input:
        mutation_data = "results/mutations/{patient}/mutation_data.csv"
    output:
        filtered = "results/analysis/{patient}/filtered_mutations.csv"
    params:
        min_coverage = config['python']['min_coverage'],
        min_freq = config['python']['min_frequency']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_mutations.py"

rule calculate_trends:
    """
    Calculate mutation trends over time.
    """
    input:
        mutations = rules.filter_mutations.output.filtered,
        metadata = "results/metadata/{patient}/metadata.csv"
    output:
        trends = "results/analysis/{patient}/mutation_trends.csv"
    params:
        p_threshold = config['python']['p_threshold']
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/calculate_trends.py"

rule map_genes:
    """
    Map mutations to genes using phylogenetic tree.
    """
    input:
        mutations = rules.filter_mutations.output.filtered,
        tree = "results/trees/{patient}/{patient}_tree.nwk"
    output:
        mapping = "results/analysis/{patient}/gene_mapping.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/map_genes.py"

rule analyze_mutations:
    """
    Analyze mutation types and patterns.
    """
    input:
        mutations = rules.filter_mutations.output.filtered,
        mapping = rules.map_genes.output.mapping
    output:
        analysis = "results/analysis/{patient}/mutation_analysis.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_mutations.py"

rule prepare_trajectories:
    """
    Prepare mutation trajectories over time.
    """
    input:
        mutations = rules.filter_mutations.output.filtered,
        metadata = "results/metadata/{patient}/metadata.csv"
    output:
        trajectories = "results/analysis/{patient}/trajectories.csv"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_trajectories.py"

rule create_plots:
    """
    Create visualization plots for analysis results.
    """
    input:
        mutations = rules.filter_mutations.output.filtered,
        trends = rules.calculate_trends.output.trends,
        mapping = rules.map_genes.output.mapping,
        analysis = rules.analyze_mutations.output.analysis,
        trajectories = rules.prepare_trajectories.output.trajectories
    output:
        plots = "results/plots/{patient}/analysis_plots.pdf"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_plots.py"

# Add all outputs to the final rule
rule all:
    input:
        filtered = expand("results/analysis/{patient}/filtered_mutations.csv", 
                         patient=config['patients']),
        trends = expand("results/analysis/{patient}/mutation_trends.csv", 
                       patient=config['patients']),
        mapping = expand("results/analysis/{patient}/gene_mapping.csv", 
                        patient=config['patients']),
        analysis = expand("results/analysis/{patient}/mutation_analysis.csv", 
                         patient=config['patients']),
        trajectories = expand("results/analysis/{patient}/trajectories.csv", 
                            patient=config['patients']),
        plots = expand("results/plots/{patient}/analysis_plots.pdf", 
                      patient=config['patients']) 