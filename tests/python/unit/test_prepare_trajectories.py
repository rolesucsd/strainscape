import pandas as pd
import tempfile
from pathlib import Path
from strainscape.prepare_trajectories import prepare_trajectories_fast

def test_prepare_trajectories_fast():
    with tempfile.TemporaryDirectory() as tmpdir:
        muts_file = Path(tmpdir) / "mut.tsv"
        trends_file = Path(tmpdir) / "trend.tsv"
        meta_file = Path(tmpdir) / "meta.csv"
        out_dir = Path(tmpdir) / "out"

        # Write minimal mutation file
        pd.DataFrame({
            'scaffold': ['chr1'],
            'position': [100],
            'ref_base': ['A'],
            'new_base': ['G'],
            'A': [10], 'C': [0], 'G': [5], 'T': [0],
            'position_coverage': [15],
            'Sample': ['S1']
        }).to_csv(muts_file, sep='\t', index=False)

        # Write minimal trends file
        pd.DataFrame({
            'scaffold': ['chr1'],
            'position': [100],
            'ref_base': ['A'],
            'new_base': ['G'],
            'slope': [0.1], 'p_value': [0.05], 'r_squared': [0.9]
        }).to_csv(trends_file, sep='\t', index=False)

        # Write minimal metadata file
        pd.DataFrame({
            'External.ID': ['S1'],
            'week_num': [1]
        }).to_csv(meta_file, index=False)

        # Run the function
        prepare_trajectories_fast(muts_file, trends_file, meta_file, out_dir)

        # Check output files exist and are not empty
        traj_file = out_dir / "mutation_trajectories.tsv"
        trends_out_file = out_dir / "mutation_trends.tsv"
        assert traj_file.exists()
        assert trends_out_file.exists()
        df = pd.read_csv(traj_file, sep='\t')
        assert not df.empty 