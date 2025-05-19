#!/usr/bin/env python3
"""
Update patients list in config.yaml based on directories in input path.
"""

import yaml
from pathlib import Path
import sys
import logging

def setup_logging():
    """Set up logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    return logging.getLogger(__name__)

def get_patient_dirs(data_dir: Path) -> list:
    """Get list of patient directories from data directory.
    
    Args:
        data_dir: Path to data directory containing patient folders
        
    Returns:
        List of patient directory names
    """
    if not data_dir.exists():
        raise FileNotFoundError(f"Data directory not found: {data_dir}")
        
    # Get all directories in data_dir
    patient_dirs = [d.name for d in data_dir.iterdir() if d.is_dir()]
    
    # Filter out any non-patient directories (e.g. logs, metadata)
    exclude_dirs = {'logs', 'metadata', '__pycache__', '.git'}
    patient_dirs = [d for d in patient_dirs if d not in exclude_dirs]
    
    return sorted(patient_dirs)

def update_config(config_path: Path, patients: list) -> None:
    """Update patients list in config.yaml.
    
    Args:
        config_path: Path to config.yaml
        patients: List of patient directory names
    """
    # Read current config
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    # Update patients list
    config['patients'] = patients
    
    # Write updated config
    with open(config_path, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

def main():
    """Main function to update patients list."""
    logger = setup_logging()
    
    # Get config path
    config_path = Path('config.yaml')
    if not config_path.exists():
        logger.error("config.yaml not found")
        sys.exit(1)
    
    # Read data directory from config
    with open(config_path) as f:
        config = yaml.safe_load(f)
    data_dir = Path(config['paths']['data_dir'])
    
    try:
        # Get patient directories
        patients = get_patient_dirs(data_dir)
        logger.info(f"Found {len(patients)} patient directories: {', '.join(patients)}")
        
        # Update config
        update_config(config_path, patients)
        logger.info(f"Updated {config_path} with patient list")
        
    except Exception as e:
        logger.error(f"Error updating patients list: {e}")
        sys.exit(1)

if __name__ == '__main__':
    main() 