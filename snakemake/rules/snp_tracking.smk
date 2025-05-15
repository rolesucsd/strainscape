"""
SNP tracking rules for the iHMP pipeline.
Handles mutation analysis and tracking over time.
"""

# Import wildcard functions
from scripts.wildcards import (
    COMBINED_SNV_INFO, COMBINED_SCAFFOLD_INFO,
    FILTERED_MUTATIONS, MUTATION_TRENDS,
    MUTATION_TRAJECTORIES, MAPPED_MUTATIONS,
    MUTATION_TYPES, BAKTA_TSV,
    ENV_INSTRAIN, ENV_PYTHON
)

# Calculate mutation trends over time
rule calculate_trends:
    input:
        mutation_file = rules.merge_snv_info.output.merged_file,
        metadata_file = config['metadata']['file']
    output:
        trends_file = os.path.join(config['paths']['results'], "mutation_trends.tsv")
    params:
        min_slope = config.get("min_slope", 0.01),
        min_abs_change = config.get("min_abs_change", 0.1),
        p_value = config.get("p_value", 0.05),
        threads = config.get("threads", 4),
        chunksize = config.get("chunksize", 10000)
    resources:
        mem_mb = 16000,
        threads = 4
    log:
        os.path.join(config['paths']['log_dir'], "calculate_trends.log")
    conda:
        config['conda_envs']['instrain']
    shell:
        """
        python {workflow.basedir}/scripts/calculate_trends.py \
            --mutation_file {input.mutation_file} \
            --metadata_file {input.metadata_file} \
            --output_file {output.trends_file} \
            --min_slope {params.min_slope} \
            --min_abs_change {params.min_abs_change} \
            --p_value {params.p_value} \
            --threads {params.threads} \
            --chunksize {params.chunksize} \
            --log_file {log} 2>&1
        """

# Prepare mutation trajectories over time
rule prepare_trajectories:
    input:
        mutations = MUTATION_TRENDS("{patient}"),
        metadata = config['metadata']['file']
    output:
        trajectories = TRAJECTORY_MUTATIONS("{patient}")
    params:
        min_trajectory_change = config["snp_tracking"]["min_trajectory_change"]
    conda:
        config['conda_envs']['instrain']
    script:
        "../scripts/prepare_trajectories.py"

# Map mutations to genes (no gene trees)
rule map_genes:
    input:
        mutations = MUTATION_TRAJECTORIES("{patient}"),
        trends = MUTATION_TRENDS("{patient}"),
        annotation = BAKTA_TSV("{patient}")
    output:
        mapped = MAPPED_MUTATIONS("{patient}")
    conda:
        config['conda_envs']['instrain']
    shell:
        """
        python {{workflow.basedir}}/scripts/map_genes.py \
            --mutation_file {{input.mutations}} \
            --trend_file {{input.trends}} \
            --gene_file {{input.annotation}} \
            --output_file {{output.mapped}}
        """

# Analyze mutation types (Silent, Missense, Nonsense)
rule analyze_mutation_types:
    input:
        mutations = MAPPED_MUTATIONS("{patient}"),
        reference = lambda wc: config["reference"][wc.patient]
    output:
        mutation_types = MUTATION_TYPES("{patient}")
    conda:
        config['conda_envs']['instrain']
    script:
        "../scripts/analyze_mutation_types.py"

# Final rule to collect all outputs
rule all:
    input:
        expand(MUTATION_TYPES("{patient}"), patient=config["patients"]) 