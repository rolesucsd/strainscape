#!/usr/bin/env python3
"""
Combine flagstat files from mapping outputs.
"""

import os
import sys
import re
from pathlib import Path

def combine_flagstats(base_dir, out_file):
    """
    Combine flagstat files from mapping outputs.
    
    Args:
        base_dir (str): Base directory containing mapping outputs
        out_file (str): Output file path for combined flagstats
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_file), exist_ok=True)
    
    # Create/truncate output file
    with open(out_file, 'w') as out:
        # Write header
        out.write("MappedFraction\tFilename\tSample\n")
        
        # Loop through each flagstat file
        base_path = Path(base_dir)
        for flagstat_file in base_path.glob("*/flagstat/*.flagstat.txt"):
            try:
                # Read first mapped line
                with open(flagstat_file, 'r') as f:
                    for line in f:
                        if 'mapped (' in line:
                            # Extract percentage
                            match = re.search(r'\(([0-9.]+)\s*%\)', line)
                            if match:
                                pct_float = float(match.group(1)) / 100.0
                                fname = flagstat_file.name
                                folder = flagstat_file.parent.parent.name
                                out.write(f"{pct_float}\t{fname}\t{folder}\n")
                            break
            except Exception as e:
                print(f"Error processing {flagstat_file}: {e}", file=sys.stderr)
                continue

def main():
    """Main function to parse arguments and run combination."""
    if len(sys.argv) != 3:
        print("Usage: combine_flagstats.py <base_dir> <out_file>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    out_file = sys.argv[2]
    
    combine_flagstats(base_dir, out_file)
    print(f"Combined flagstat metrics written to: {out_file}")

if __name__ == "__main__":
    main() 