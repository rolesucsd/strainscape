"""
SNP tracking rules for the iHMP pipeline.
Handles mutation analysis and tracking over time.
"""
rule process_scaffolds:
    """
    Process combined scaffold information for downstream analysis.
    
    This rule processes the combined scaffold information by:
    1. Filtering scaffolds based on length, coverage, and breadth thresholds
    2. Merging with STB and metadata information
    3. Calculating summary statistics
    
    Input:
        scaffold_file: Combined scaffold information file
        bin_file: Path to bin.txt file (scaffold-bin mapping)
    
    Output:
        processed_file: Processed scaffold information file
    
    Parameters:
        min_length: Minimum scaffold length (default: 1000)
        min_coverage: Minimum coverage threshold (default: 10.0)
        min_breadth: Minimum breadth threshold (default: 0.4)
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        scaffold_file = COMBINED_SCAFFOLD_INFO("{patient}"),
        bin_file = BIN_TXT("{patient}")
    output:
        processed_scaffolds = PROCESSED_SCAFFOLD("{patient}")
    params:
        min_length = 1000,
        min_coverage = 5,
        min_breadth = 0.4
    shell:
        """
        python strainscape/process_scaffolds.py \
            --scaffold_file {input.scaffold_file} \
            --bin_file {input.bin_file} \
            --output_file {output.processed_scaffolds} \
            --min_length {params.min_length} \
            --min_coverage {params.min_coverage} \
            --min_breadth {params.min_breadth}
        """

rule merge_snv_info:
    """
    Merge SNV information with metadata and filtered scaffolds.
    
    This rule combines SNV information with metadata and filtered scaffold information
    to ensure SNVs only map to scaffolds that passed previous filtering criteria.
    
    Input:
        snv_file: Combined SNV information file
        metadata_file: Metadata file with sample information
        processed_scaffold_file: Processed scaffold file containing filtered scaffolds
    
    Output:
        filtered_file: Filtered SNV information file
    
    Parameters:
        min_coverage: Minimum coverage threshold (default: 10)
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        snv_file = COMBINED_SNV_INFO("{patient}"),
        processed_scaffolds_file = PROCESSED_SCAFFOLD("{patient}")
    output:
        filtered_file = FILTERED_MUTATIONS("{patient}")
    params:
        metadata_file = config['metadata']['file'],
        min_coverage = 10
    shell:
        """
        python strainscape/filter_mutations.py \
            --snv-file {input.snv_file} \
            --output-file {output.filtered_file} \
            --metadata-file {params.metadata_file} \
            --processed-scaffolds-file {input.processed_scaffolds_file} \
            --min-coverage {params.min_coverage} \
            --log-file {log} 2>&1
        """

rule calculate_trends:
    """
    Calculate mutation trends and trajectories over time.
    
    This rule analyzes mutation frequencies across timepoints to identify
    significant trends. It uses the following parameters:
    - Minimum slope: 0.01 (default)
    - P-value threshold: 0.05 (default)
    
    Input:
        filtered_file: File containing filtered mutations
        metadata_file: Metadata file with sample information
    
    Output:
        trends_file: File containing significant mutation trends
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        filtered_file = FILTERED_MUTATIONS("{patient}")
    output:
        trends_file = MUTATION_TRENDS("{patient}")
    params:
        min_slope = 0.01,
        p_value = 0.05
    shell:
        """
        python strainscape/calculate_trends.py \
            --mutation_file {input.filtered_file} \
            --output_file {output.trends_file} \
            --min_slope {params.min_slope} \
            --p_value {params.p_value}
        """

rule map_genes:
    """
    Map mutations to genes using BAKTA annotations.
    
    This rule maps mutations to genes using the BAKTA gene annotations.
    It identifies which mutations fall within coding regions and
    determines their potential impact.
    
    Input:
        trajectory_file: File containing mutation trajectories
        gene_file: BAKTA TSV file with gene annotations
    
    Output:
        mapped_file: File containing mutations mapped to genes
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        trends_file = MUTATION_TRENDS("{patient}"),
        gene_file = BAKTA_TSV("{patient}")
    output:
        mapped_file = MAPPED_MUTATIONS("{patient}")
    shell:
        """
        python strainscape/map_genes.py \
            --trend_file {input.trends_file} \
            --gene_file {input.gene_file}
        """

rule analyze_mutation_types:
    """
    Analyze types of mutations in mapped genes.
    
    This rule analyzes the types of mutations (e.g., synonymous, non-synonymous)
    in the mapped genes using the BAKTA gene sequences.
    
    Input:
        fa_file: BAKTA FNA file with gene sequences
        mapped_file: File containing mapped mutations
    
    Output:
        analyzed_file: File containing mutation type analysis
    
    Resources:
        mem_mb: 16000
        threads: 4
    """
    input:
        fa_file = BAKTA_FNA("{patient}"),
        mapped_file = MAPPED_MUTATIONS("{patient}")
    output:
        analyzed_file = MUTATION_TYPES("{patient}")
    shell:
        """
        python strainscape/analyze_mutation_types.py \
            --mutation_file {input.mapped_file} \
            --reference_file {input.fa_file} \
            --output_file {output.analyzed_file}
        """ 