"""
Rules for creating phylogenetic trees from Bakta annotation files.
"""

rule create_trees:
    """
    Create phylogenetic trees from Bakta annotation files.
    """
    input:
        bakta_tsv = "results/annotation/{sample}/bakta/{sample}.tsv"
    output:
        tree = "results/trees/{sample}/{sample}_tree.nwk",
        tree_stats = "results/trees/{sample}/{sample}_tree_stats.txt"
    params:
        min_coverage = config['python']['min_coverage'],
        min_quality = config['python']['min_quality']
    conda:
        "../envs/python.yaml"
    shell:
        """
        python src/python/ihmp/analysis/create_trees.py \
            --input {input.bakta_tsv} \
            --output {output.tree} \
            --stats {output.tree_stats} \
            --min-coverage {params.min_coverage} \
            --min-quality {params.min_quality}
        """

# Add tree outputs to the final rule
rule all:
    input:
        trees = expand("results/trees/{sample}/{sample}_tree.nwk", 
                      sample=config['samples']),
        tree_stats = expand("results/trees/{sample}/{sample}_tree_stats.txt", 
                          sample=config['samples']) 