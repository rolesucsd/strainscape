#!/usr/bin/env python3
"""
Utility functions for the strainscape package.

This module provides common utility functions used across the strainscape package,
including logging setup, file operations, and data processing helpers.
"""

import logging
import os
from pathlib import Path
from typing import Dict, List, Optional, Union
import yaml
import pandas as pd
from functools import lru_cache
import time

def setup_logging(log_file: Optional[Path] = None) -> None:
    """Set up logging configuration.
    
    Args:
        log_file: Optional path to write log file. If None, logs to stderr.
        
    Returns:
        None. Configures the root logger with appropriate handlers.
    """
    # Get the root logger
    logger = logging.getLogger()
    
    # If logger already has handlers, don't add more
    if logger.handlers:
        return
        
    # Set basic configuration
    logger.setLevel(logging.INFO)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Add console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # Add file handler if log_file specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for a module.
    
    Args:
        name: Name of the module (typically __name__)
        
    Returns:
        Logger instance configured with the module's name
    """
    return logging.getLogger(name)

def ensure_dir(directory: Path) -> None:
    """Ensure a directory exists, creating it if necessary.
    
    Args:
        directory: Path to the directory to ensure exists
        
    Returns:
        None. Creates the directory if it doesn't exist.
    """
    directory.mkdir(parents=True, exist_ok=True)

def safe_divide(a: float, b: float) -> float:
    """Safely divide two numbers, returning 0 if denominator is 0.
    
    Args:
        a: Numerator
        b: Denominator
        
    Returns:
        a/b if b != 0, else 0
    """
    return a / b if b != 0 else 0

def format_percentage(value: float) -> str:
    """Format a float as a percentage string.
    
    Args:
        value: Float to format (e.g. 0.123 for 12.3%)
        
    Returns:
        String formatted as percentage with 1 decimal place
    """
    return f"{value * 100:.1f}%"

def parse_float(value: str) -> Optional[float]:
    """Parse a string to float, returning None if invalid.
    
    Args:
        value: String to parse
        
    Returns:
        Float value if valid, None if invalid
    """
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

def parse_int(value: str) -> Optional[int]:
    """Parse a string to integer, returning None if invalid.
    
    Args:
        value: String to parse
        
    Returns:
        Integer value if valid, None if invalid
    """
    try:
        return int(value)
    except (ValueError, TypeError):
        return None

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

# Data loading and saving
def load_data(file_path: str) -> pd.DataFrame:
    """Load data from a CSV file.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        DataFrame containing the loaded data
    """
    return pd.read_csv(file_path, dtype=str)

def save_data(data: pd.DataFrame, file_path: str) -> None:
    """Save data to a CSV file.
    
    Args:
        data: DataFrame to save
        file_path: Path where to save the file
    """
    # Ensure the directory exists
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    data.to_csv(file_path, index=False)

def filter_data(data: pd.DataFrame, condition: str) -> pd.DataFrame:
    """
    Filter data based on a condition.
    
    Args:
        data: DataFrame to filter
        condition: String containing the filter condition
    
    Returns:
        Filtered DataFrame
    """
    return data.query(condition)

def merge_data(df1: pd.DataFrame, df2: pd.DataFrame, on: str, how: str = 'inner') -> pd.DataFrame:
    """Merge two DataFrames.
    
    Args:
        df1: First DataFrame
        df2: Second DataFrame
        on: Column to merge on
        how: Type of merge to perform
        
    Returns:
        Merged DataFrame
    """
    return pd.merge(df1, df2, on=on, how=how)

def calculate_statistics(data: pd.DataFrame, column: str, groupby: Optional[str] = None) -> Union[Dict, pd.DataFrame]:
    """
    Calculate basic statistics for a column.
    
    Args:
        data: DataFrame containing the data
        column: Column to calculate statistics for
        groupby: Optional column to group by
    
    Returns:
        Dictionary of statistics or DataFrame with grouped statistics
    """
    if groupby:
        return data.groupby(groupby)[column].agg(['mean', 'std', 'min', 'max']).reset_index()
    else:
        return {
            'mean': data[column].mean(),
            'std': data[column].std(),
            'min': data[column].min(),
            'max': data[column].max()
        }

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

# Create root logger instance
logger = setup_logging()

def get_output_path(base_path, patient, sample=None, suffix=None):
    """Generate output path for a file."""
    path = Path(base_path) / patient
    if sample:
        path = path / sample
    if suffix:
        path = path.with_suffix(suffix)
    return str(path) 