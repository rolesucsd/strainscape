#!/usr/bin/env python3
"""
Validate input files and dependencies for the iHMP pipeline.
This script checks for the existence and validity of all required input files
and configurations before running the pipeline.
"""

import os
import sys
import yaml
import logging
from pathlib import Path

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def load_config(config_file):
    """Load and validate configuration file."""
    try:
        with open(config_file, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        raise ValueError(f"Error loading config file: {e}")

def check_directory(path, name):
    """Check if directory exists and is accessible."""
    if not os.path.exists(path):
        raise ValueError(f"{name} directory not found: {path}")
    if not os.access(path, os.R_OK):
        raise ValueError(f"{name} directory not readable: {path}")

def check_file(path, name):
    """Check if file exists and is accessible."""
    if not os.path.exists(path):
        raise ValueError(f"{name} file not found: {path}")
    if not os.access(path, os.R_OK):
        raise ValueError(f"{name} file not readable: {path}")

def validate_reference_data(config, logger):
    """Validate reference data files and directories."""
    ref = config.get('reference', {})
    required_ref = ['bakta_db', 'genomes', 'annotations', 'taxonomy']
    
    for item in required_ref:
        path = ref.get(item)
        if not path:
            raise ValueError(f"Missing reference {item} path in config")
        if item == 'taxonomy':
            check_file(path, f"Reference {item}")
        else:
            check_directory(path, f"Reference {item}")

def validate_metadata(config, logger):
    """Validate metadata files."""
    metadata = config.get('metadata', {})
    required_meta = ['file', 'patients']
    
    for item in required_meta:
        path = metadata.get(item)
        if not path:
            raise ValueError(f"Missing metadata {item} path in config")
        check_file(path, f"Metadata {item}")

def validate_samples(config, logger):
    """Validate sample files and directories."""
    samples_file = os.path.join(os.path.dirname(config['samples']), 'samples.yaml')
    patients_file = os.path.join(os.path.dirname(config['patients']), 'patients.yaml')
    
    check_file(samples_file, "Samples")
    check_file(patients_file, "Patients")
    
    # Load and validate sample data
    with open(samples_file, 'r') as f:
        samples = yaml.safe_load(f)
    
    for patient, samples in samples.items():
        for sample, files in samples.items():
            for file_type, file_path in files.items():
                check_file(file_path, f"Sample file {file_type} for {patient}/{sample}")

def validate_gene_trees(config, logger):
    """Validate gene trees configuration."""
    gene_trees = config.get('gene_trees', {})
    if not gene_trees.get('directory'):
        raise ValueError("Missing gene trees directory in config")
    check_directory(gene_trees['directory'], "Gene trees")

def main():
    """Main function to validate all inputs."""
    logger = setup_logging()
    
    if len(sys.argv) != 2:
        logger.error("Usage: validate_inputs.py <config_file>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    logger.info(f"Loading configuration from {config_file}")
    
    try:
        config = load_config(config_file)
        
        # Validate all components
        logger.info("Validating reference data...")
        validate_reference_data(config, logger)
        
        logger.info("Validating metadata...")
        validate_metadata(config, logger)
        
        logger.info("Validating samples...")
        validate_samples(config, logger)
        
        logger.info("Validating gene trees...")
        validate_gene_trees(config, logger)
        
        logger.info("All validations passed successfully!")
        
    except Exception as e:
        logger.error(f"Validation failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main() 