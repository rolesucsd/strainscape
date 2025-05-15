"""
Utility functions and logging configuration for the iHMP pipeline.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import yaml
import pandas as pd
from functools import lru_cache
import time

# Configure logging
def setup_logging(log_file: Optional[str] = None, level: int = logging.INFO) -> logging.Logger:
    """
    Set up logging configuration for the pipeline.
    
    Args:
        log_file: Optional path to log file. If None, logs to console only.
        level: Logging level (default: INFO)
    
    Returns:
        Logger instance
    """
    logger = logging.getLogger('ihmp_pipeline')
    logger.setLevel(level)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if log_file specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

# Performance monitoring
class PerformanceMonitor:
    """Monitor and log performance metrics."""
    
    def __init__(self, logger: logging.Logger):
        self.logger = logger
        self.metrics: Dict[str, float] = {}
    
    def start_operation(self, operation: str) -> None:
        """Start timing an operation."""
        self.metrics[operation] = time.time()
    
    def end_operation(self, operation: str) -> None:
        """End timing an operation and log duration."""
        if operation in self.metrics:
            duration = time.time() - self.metrics[operation]
            self.logger.info(f"Operation '{operation}' took {duration:.2f} seconds")
            del self.metrics[operation]

# Data validation
def validate_file_exists(file_path: Union[str, Path]) -> bool:
    """
    Validate that a file exists and is readable.
    
    Args:
        file_path: Path to the file to validate
    
    Returns:
        bool: True if file exists and is readable
    
    Raises:
        FileNotFoundError: If file doesn't exist
        PermissionError: If file isn't readable
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"File not readable: {file_path}")
    return True

@lru_cache(maxsize=128)
def load_config(config_file: str) -> Dict:
    """
    Load and cache configuration file.
    
    Args:
        config_file: Path to config file
    
    Returns:
        Dict containing configuration
    """
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def get_samples(patient: str) -> Dict[str, str]:
    """
    Get sample information for a patient.
    
    Args:
        patient: Patient identifier
    
    Returns:
        Dict mapping sample names to file paths
    """
    config = load_config("config/config.yaml")
    return config["samples"][patient]

# Memory efficient file reading
def read_large_csv(file_path: str, chunksize: int = 10000) -> pd.DataFrame:
    """
    Read large CSV files efficiently using chunks.
    
    Args:
        file_path: Path to CSV file
        chunksize: Number of rows to read at once
    
    Returns:
        DataFrame containing all data
    """
    chunks = []
    for chunk in pd.read_csv(file_path, chunksize=chunksize):
        chunks.append(chunk)
    return pd.concat(chunks, ignore_index=True)

# Create logger instance
logger = setup_logging()

def get_output_path(base_path, patient, sample=None, suffix=None):
    """Generate output path for a file."""
    path = Path(base_path) / patient
    if sample:
        path = path / sample
    if suffix:
        path = path.with_suffix(suffix)
    return str(path) 