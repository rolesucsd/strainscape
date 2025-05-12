import glob
import os
import re
import pandas as pd

def build_mega_mutation_df(root_dir="/ddn_scratch/roles/strain_analysis/iHMP", 
                           pattern="SNV_filtered_trend_trajectory_mapped_mutation.txt"):
    """
    Walk all patient subdirs under root_dir, read each SNV file, filter on OLS_pvalue ≤ 0.05,
    pivot numeric columns to long format, and concatenate everything into one DataFrame.
    """
    # 1) find all files
    file_paths = glob.glob(os.path.join(root_dir, "*", pattern))
    all_long = []

    for fp in file_paths:
        # 2) derive sample name from the patient directory
        sample = os.path.basename(os.path.dirname(fp))

        # 3) read into pandas (everything as string so we can pivot safely)
        df = pd.read_csv(fp, sep="\t", dtype=str, na_values=["", "NA"])
        df["Sample"] = sample

        # 4) filter on OLS_pvalue ≤ 0.1(coerce to numeric for comparison)
        df["OLS_pvalue_num"] = pd.to_numeric(df["OLS_pvalue"], errors="coerce")
        df = df[df["OLS_pvalue_num"] <= 0.1]
        df.drop(columns="OLS_pvalue_num", inplace=True)

        # 5) identify “numeric” columns by name (all digits or digits+dot)
        numeric_cols = [c for c in df.columns 
                        if re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", c)]
        # ensure they’re strings (to mimic R’s as.character)
        df[numeric_cols] = df[numeric_cols].astype(str)

        # 6) pivot longer
        id_vars = [c for c in df.columns if c not in numeric_cols]
        long_df = df.melt(id_vars   = id_vars,
                          value_vars= numeric_cols,
                          var_name   = "week",
                          value_name = "percentage")

        print(f" Pivoted sample {sample}: {len(numeric_cols)} weeks → {long_df.shape[0]} rows")
        all_long.append(long_df)

    # 7) concatenate (pandas will union the columns, filling missing with NaN)
    mega = pd.concat(all_long, ignore_index=True, sort=False)
    # (optional) cast all columns back to string
    mega = mega.astype(str, errors="ignore")

    print(f"\nCombined into mega DataFrame: {mega.shape[0]} rows × {mega.shape[1]} columns")
    return mega

# Example usage in your main script:
if __name__ == "__main__":
    mega_mutations = build_mega_mutation_df()
    # write out if you like:
    mega_mutations.to_csv("all_patients_mutations_mega.tsv", 
                          sep="\t", index=False)
