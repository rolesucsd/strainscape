"""
Utility functions for the iHMP pipeline.
"""

import os
import yaml
from pathlib import Path

def load_config(config_file):
    """Load configuration from YAML file."""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def get_samples(patient):
    """Get sample paths for a patient."""
    config = load_config("config/config.yaml")
    raw_data = Path(config["paths"]["raw_data"])
    patient_dir = raw_data / patient
    
    if not patient_dir.exists():
        raise ValueError(f"Patient directory not found: {patient_dir}")
    
    samples = {}
    for sample_dir in patient_dir.iterdir():
        if sample_dir.is_dir():
            r1 = list(sample_dir.glob("*_R1_*.fastq.gz"))
            if r1:
                samples[sample_dir.name] = str(r1[0])
    
    return samples

def get_output_path(base_path, patient, sample=None, suffix=None):
    """Generate output path for a file."""
    path = Path(base_path) / patient
    if sample:
        path = path / sample
    if suffix:
        path = path.with_suffix(suffix)
    return str(path) 