#!/usr/bin/env python3
# calculate_combined_stats_v9.py --------------------------------------------
# v8 + convergent-evolution analysis  (KO enrichment on sweeping variants)
# ---------------------------------------------------------------------------

import argparse, re, sys, gzip, uuid, functools as ft, random, math
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds
from tqdm import tqdm

pd.options.mode.copy_on_write = True           # pandas ≥2.1

# ───────────────────────────── CLI ──────────────────────────────────────────
def get_args():
    p = argparse.ArgumentParser("iHMP mutation summariser (stream, sweeps, KO-convergence)")
    p.add_argument("--mut-file",      required=True)
    p.add_argument("--scaffold-file", required=True)
    p.add_argument("--meta-file",     required=True)
    p.add_argument("--out-dir",       default="summaries")
    p.add_argument("--chunk-size",    type=int, default=1_000_000)
    p.add_argument("--low-mem",       action="store_true")
    p.add_argument("--patient-chunks",action="store_true")
    # sweep windows
    p.add_argument("--base-window",   default="-2:0")
    p.add_argument("--post-window",   default="10:28")
    # convergent-evolution
    p.add_argument("--max-lg-size",   type=int, default=100,
                   help="analyse LGs ≤ this many variants")
    p.add_argument("--permutations",  type=int, default=100_000,
                   help="number of permutations for KO empirical P values")
    p.add_argument("--lg-col",        default="lg",
                   help="column that contains linkage-group IDs")
    return p.parse_args()

def parse_win(win:str)->Tuple[int,int]:
    a,b = map(int,win.split(":")); return (min(a,b), max(a,b))

# ────────────────────────── helpers ─────────────────────────────────────────
detect_week_cols = lambda cols: sorted(int(c) for c in cols if re.fullmatch(r"\d+",c))
def classify_dn_ds(mut_type):
    m = mut_type.str.lower().fillna("")
    return m.isin(["missense","nonsense"]), m.eq("silent")
split_dbxrefs=lambda s:[t.strip().lower() for t in str(s).split(",") if t.strip()]
extract_go   =lambda tok:[t for t in tok if t.startswith("go:")]
extract_kegg =lambda tok:[t for t in tok if t.startswith("ec:") or re.fullmatch(r"k\d{5}",t)]

# ───────────────────────── bin map ──────────────────────────────────────────
def load_bin_map(scaffold_file, meta_file):
    meta = pd.read_csv(meta_file, low_memory=False)
    meta.columns = meta.columns.str.strip().str.lower().str.replace(" ","_")
    sample_col = "external.id" if "external.id" in meta.columns else "external_id"
    meta = meta.rename(columns={sample_col:"sample"})[["sample","participant_id"]]\
             .rename(columns={"participant_id":"patient_id"}).dropna()

    scaf = pd.read_csv(scaffold_file, sep="\t",
                       usecols=["scaffold","Sample","bin"], low_memory=False)
    scaf.columns = scaf.columns.str.lower()
    mp = scaf.merge(meta, on="sample", how="left").dropna(subset=["patient_id"])
    return {(r.patient_id, r.scaffold): r.bin for r in mp.itertuples(index=False)}

# ───────────────────── streaming writer ─────────────────────────────────────
def append_parquet(dir_path:Path, df):
    if df.empty: return
    dir_path.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df, preserve_index=False),
                   dir_path/f"part-{uuid.uuid4().hex}.parquet",
                   compression="snappy")

# ─────────── incremental aggregation (counts & sums) ───────────────────────
def merge_stat(running,new,key_cols,agg_dict):
    if new.empty: return running
    combined = pd.concat([running,new]) if running is not None else new
    return combined.groupby(key_cols,as_index=False).agg(agg_dict)

# ─────────────────────── sweep helpers ─────────────────────────────────────
def tag_window(day, base_rng, post_rng, is_control=False):
    lo1,hi1 = (0,2)   if is_control else base_rng
    lo2,hi2 = (49,63) if is_control else post_rng
    if lo1<=day<=hi1: return "baseline"
    if lo2<=day<=hi2: return "post"
    return None

# ─────────────────────── chunk processor ────────────────────────────────────
def process_chunk(df, week_cols, bin_map, base_rng, post_rng,
                  control_ids:set, lg_col:str):
    # numeric coercion
    num_cols = ["p_value","slope","freq_range","min_freq","max_freq","r_sq",
                *map(str,week_cols)]
    for c in num_cols:
        if c in df.columns: df[c] = pd.to_numeric(df[c],errors="coerce")

    # bin
    if "bin" not in df.columns:
        df["bin"]=[bin_map.get((pid,chrom)) for pid,chrom in zip(df.patient_id,df.chromosome)]

    # DbXrefs → GO / KEGG
    if "dbxrefs" in df.columns and "go_terms" not in df.columns:
        toks = df.dbxrefs.apply(split_dbxrefs)
        df["go_terms"]=toks.apply(extract_go)
        df["kegg_terms"]=toks.apply(extract_kegg)

    # ───────────── long freq (burden + sweeps) ─────────────
    id_cols=[c for c in df.columns if c not in map(str,week_cols)]
    long_df=(df.melt(id_vars=id_cols,value_vars=list(map(str,week_cols)),
                     var_name="week_num",value_name="freq").dropna(subset=["freq"]))
    long_df["week_num"]=long_df.week_num.astype(int)
    long_df["window"]=[tag_window(w,base_rng,post_rng,pid in control_ids)
                       for pid,w in zip(long_df.patient_id,long_df.week_num)]

    # sweep sums
    sweep_df=(long_df[long_df.window.notna()]
                .groupby(["patient_id","chromosome","position","window"],
                         as_index=False)
                .agg(sum_freq=("freq","sum"),n=("freq","size")))

    # variant meta (for convergent evolution) – one row per variant
    meta_cols=["patient_id","chromosome","position","gene","kegg_terms",lg_col,"bin"]
    var_meta=df.drop_duplicates(subset=["patient_id","chromosome","position"])[meta_cols]

    # per-SNV stats
    snp_cols=[c for c in ["patient_id","chromosome","position","ref_base","new_base",
                          "slope","p_value","r_sq","min_freq","max_freq","freq_range"]
              if c in df.columns]
    snp_stats=df[snp_cols].drop_duplicates()

    # gene trajectory (sums not means)
    gene_traj=pd.DataFrame()
    if {"gene","p_value"}.issubset(df.columns):
        gene_traj=(df.groupby(["patient_id","gene"],as_index=False)
                     .agg(n_snps=("p_value","size"),
                          min_p=("p_value","min"),
                          sum_p=("p_value","sum"),
                          sum_abs_slope=("slope",lambda x:np.nansum(np.abs(x))),
                          sum_delta=("freq_range","sum")))

    # GO / KEGG
    go_stats=kegg_stats=pd.DataFrame(),pd.DataFrame()
    if "go_terms" in df.columns:
        go_exp=(df[["patient_id","go_terms","p_value","slope","freq_range"]]
                  .explode("go_terms").dropna(subset=["go_terms"]))
        if not go_exp.empty:
            go_stats=(go_exp.groupby(["patient_id","go_terms"],as_index=False)
                        .agg(n_snps=("p_value","size"),min_p=("p_value","min"),
                             sum_p=("p_value","sum"),
                             sum_abs_slope=("slope",lambda x:np.nansum(np.abs(x))),
                             sum_delta=("freq_range","sum")))
        kegg_exp=(df[["patient_id","kegg_terms","p_value","slope","freq_range"]]
                    .explode("kegg_terms").dropna(subset=["kegg_terms"]))
        if not kegg_exp.empty:
            kegg_stats=(kegg_exp.groupby(["patient_id","kegg_terms"],as_index=False)
                          .agg(n_snps=("p_value","size"),min_p=("p_value","min"),
                               sum_p=("p_value","sum"),
                               sum_abs_slope=("slope",lambda x:np.nansum(np.abs(x))),
                               sum_delta=("freq_range","sum")))

    # dN/dS
    gene_dnds=bin_dnds=pd.DataFrame(),pd.DataFrame()
    if "mutation_type" in df.columns:
        is_dn,is_ds=classify_dn_ds(df.mutation_type)
        if "gene" in df.columns:
            gene_dnds=(df.assign(is_dn=is_dn,is_ds=is_ds)
                         .groupby(["patient_id","gene"],as_index=False)
                         .agg(dn=("is_dn","sum"),ds=("is_ds","sum")))
        if "bin" in df.columns:
            bin_dnds=(df.assign(is_dn=is_dn,is_ds=is_ds)
                        .dropna(subset=["bin"])
                        .groupby(["patient_id","bin"],as_index=False)
                        .agg(dn=("is_dn","sum"),ds=("is_ds","sum")))

    return (long_df,snp_stats,gene_traj,gene_dnds,bin_dnds,
            go_stats,kegg_stats,sweep_df,var_meta)

# ─────────────────────────── main ───────────────────────────────────────────
def main():
    a=get_args(); out_dir=Path(a.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    base_rng=parse_win(a.base_window); post_rng=parse_win(a.post_window)

    meta=pd.read_csv(a.meta_file,low_memory=False)
    meta.columns=meta.columns.str.strip().str.lower().str.replace(" ","_")
    control_ids=set(meta.loc[meta.antibiotics.isna(),"participant_id"])

    bin_map=load_bin_map(a.scaffold_file,a.meta_file)

    preview=pd.read_csv(a.mut_file,sep="\t",nrows=2000,low_memory=False)
    week_cols=detect_week_cols(preview.columns)
    if not week_cols: sys.exit("❌ No numeric week columns detected.")
    print(f"Detected {len(week_cols)} week columns; first 10 → {week_cols[:10]}")

    opener=gzip.open if a.mut_file.endswith(".gz") else open
    total_rows=sum(1 for _ in opener(a.mut_file,"rt"))-1
    pbar=tqdm(total=total_rows,unit="rows")

    runs={k:None for k in ["gene_traj","gene_dnds","bin_dnds",
                           "go_stats","kegg_stats","sweep_core"]}
    # we stream-write three heavy dirs
    dirs={"mutation_long","snp_stats","variant_meta"}
    for d in dirs: (out_dir/d).mkdir(parents=True, exist_ok=True)

    for chunk in pd.read_csv(a.mut_file, sep="\t", chunksize=a.chunk_size,
                             low_memory=not a.low_mem):
        chunk.columns=(chunk.columns.str.strip().str.lower().str.replace(" ","_"))
        sub_chunks=(chunk.groupby("participant_id",sort=False)
                          .apply(lambda g:g)) if a.patient_chunks else [chunk]
        for sub in sub_chunks:
            if sub.empty: continue
            if "participant_id" in sub.columns:
                sub=sub.rename(columns={"participant_id":"patient_id"})
            if sub.patient_id.isna().all(): continue

            dfs=process_chunk(sub,week_cols,bin_map,base_rng,post_rng,
                              control_ids,a.lg_col)
            (long_df,snp_stats,gene_traj,gene_dnds,bin_dnds,
             go_stats,kegg_stats,sweep_df,var_meta)=dfs

            append_parquet(out_dir/"mutation_long", long_df)
            append_parquet(out_dir/"snp_stats",     snp_stats)
            append_parquet(out_dir/"variant_meta",  var_meta)

            runs["sweep_core"]=merge_stat(
                runs["sweep_core"], sweep_df,
                ["patient_id","chromosome","position","window"],
                {"sum_freq":"sum","n":"sum"})

            runs["gene_traj"]=merge_stat(runs["gene_traj"],gene_traj,
                ["patient_id","gene"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["gene_dnds"]=merge_stat(runs["gene_dnds"],gene_dnds,
                ["patient_id","gene"],{"dn":"sum","ds":"sum"})
            runs["bin_dnds"]=merge_stat(runs["bin_dnds"],bin_dnds,
                ["patient_id","bin"],{"dn":"sum","ds":"sum"})
            runs["go_stats"]=merge_stat(runs["go_stats"],go_stats,
                ["patient_id","go_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
            runs["kegg_stats"]=merge_stat(runs["kegg_stats"],kegg_stats,
                ["patient_id","kegg_terms"],
                {"n_snps":"sum","min_p":"min",
                 "sum_p":"sum","sum_abs_slope":"sum","sum_delta":"sum"})
        pbar.update(len(chunk))
    pbar.close()

    # ── finalise sweep stats ────────────────────────────────────────────────
    core=runs["sweep_core"]
    sweep_stats=None
    if core is not None:
        pivot=(core.pivot(index=["patient_id","chromosome","position"],
                          columns="window",values=["sum_freq","n"]).fillna(0))
        pivot.columns=["_".join(c).strip() for c in pivot.columns]
        pivot["baseline_f"]=pivot["sum_freq_baseline"]/pivot["n_baseline"].replace(0,np.nan)
        pivot["post_f"]    =pivot["sum_freq_post"]    /pivot["n_post"].replace(0,np.nan)
        pivot[["baseline_f","post_f"]]=pivot[["baseline_f","post_f"]].fillna(0)
        rev=pivot["baseline_f"]>.5
        pivot.loc[rev,["baseline_f","post_f"]]=1-pivot.loc[rev,["baseline_f","post_f"]]
        pivot["sweeps"]=(pivot["baseline_f"]<.20)&(pivot["post_f"]>.80)
        sweep_stats=pivot.reset_index()
        pq.write_table(pa.Table.from_pandas(sweep_stats,preserve_index=False),
                       out_dir/"sweep_stats.parquet",compression="snappy")
        print(f"✔ sweep_stats.parquet ({sweep_stats.shape[0]:,} variants)")

    # ── convergent evolution  (KO enrichment) ───────────────────────────────
    if sweep_stats is not None and (out_dir/"variant_meta").exists():
        # bring in meta
        varmeta_ds=ds.dataset(out_dir/"variant_meta",format="parquet")
        meta_df=varmeta_ds.to_table().to_pandas()
        sweep_df=sweep_stats.merge(meta_df,
                                   on=["patient_id","chromosome","position"],
                                   how="left")
        # keep sweeping variants w/ LG present
        lgcol=a.lg_col
        sweep_df=sweep_df[(sweep_df.sweeps) & (~sweep_df[lgcol].isna())]

        # variant count per LG instance
        lg_counts=(sweep_df.groupby(["patient_id",lgcol],as_index=False)
                             .size().rename(columns={"size":"n_var"}))
        small_lg=set(tuple(r) for r in lg_counts.loc[lg_counts.n_var<=a.max_lg_size,
                                                     ["patient_id",lgcol]].itertuples(index=False))
        mask=sweep_df.apply(lambda r:(r.patient_id,r[lgcol]) in small_lg,axis=1)
        sweep_df=sweep_df[mask]

        # weight per LG instance =1
        ko_weights={}
        ko_participants={}
        for (pid,lg), grp in sweep_df.groupby(["patient_id",lgcol]):
            genes=grp["gene"].dropna().unique()
            if len(genes)==0: continue
            gene_w=1/len(genes)
            for g in genes:
                kos=set(ft.reduce(op.add,
                                  grp.loc[grp.gene==g,"kegg_terms"]
                                         .dropna()
                                         .apply(lambda x:x if isinstance(x,list) else [x])
                                         .tolist(),
                                  []))
                if not kos: continue
                ko_w=gene_w/len(kos)
                for k in kos:
                    ko_weights[k]=ko_weights.get(k,0)+ko_w
                    ko_participants.setdefault(k,set()).add(pid)

        obs=pd.DataFrame({"ko":list(ko_weights.keys()),
                          "observed_w":list(ko_weights.values())})
        obs["n_participants"]=obs.ko.map(lambda k: len(ko_participants[k]))

        # --- permutation background ---------------------------------------
        # pre-assemble gene universe per (patient, LG) bin (=MAG)
        bin_gene={}
        for (pid,b), grp in sweep_df.groupby(["patient_id","bin"]):
            bin_gene.setdefault((pid,b), set()).update(grp["gene"].dropna().unique())

        rng=np.random.default_rng(seed=1)
        perm_tot={k:0.0 for k in ko_weights}
        P=a.permutations
        for i in tqdm(range(P),desc="Permuting KO weights",unit="perm"):
            perm_w={k:0.0 for k in ko_weights}
            for (pid,lg), grp in sweep_df.groupby(["patient_id",lgcol]):
                b=grp["bin"].iloc[0]
                genes=list(grp["gene"].dropna().unique())
                if not genes: continue
                pool=list(bin_gene.get((pid,b),[]))
                if len(pool)<len(genes): continue
                rnd_genes=rng.choice(pool,size=len(genes),replace=False)
                gene_w=1/len(rnd_genes)
                for g in rnd_genes:
                    kos=set(ft.reduce(op.add,
                                      grp.loc[grp.gene==g,"kegg_terms"]
                                         .dropna()
                                         .apply(lambda x:x if isinstance(x,list) else [x])
                                         .tolist(),
                                      []))
                    if not kos: continue
                    kw=gene_w/len(kos)
                    for k in kos:
                        perm_w[k]+=kw
            for k in perm_w: perm_tot[k]+=perm_w[k]
        mean_exp={k:perm_tot[k]/P for k in ko_weights}

        # empirical P value
        geq_counts={k:0 for k in ko_weights}
        # second pass to avoid storing all perms
        rng=np.random.default_rng(seed=2)
        for i in tqdm(range(P),desc="Empirical P",unit="perm"):
            perm_w={k:0.0 for k in ko_weights}
            for (pid,lg), grp in sweep_df.groupby(["patient_id",lgcol]):
                b=grp["bin"].iloc[0]
                genes=list(grp["gene"].dropna().unique())
                if not genes: continue
                pool=list(bin_gene.get((pid,b),[]))
                if len(pool)<len(genes): continue
                rnd_genes=rng.choice(pool,size=len(genes),replace=False)
                gene_w=1/len(rnd_genes)
                for g in rnd_genes:
                    kos=set(ft.reduce(op.add,
                                      grp.loc[grp.gene==g,"kegg_terms"]
                                         .dropna()
                                         .apply(lambda x:x if isinstance(x,list) else [x])
                                         .tolist(),
                                      []))
                    if not kos: continue
                    kw=gene_w/len(kos)
                    for k in kos: perm_w[k]+=kw
            for k in geq_counts:
                if perm_w[k]>=ko_weights[k]-1e-12: geq_counts[k]+=1
        obs["mean_exp_w"]=obs.ko.map(mean_exp)
        obs["p_emp"]=obs.ko.map(lambda k:(geq_counts[k]+1)/(P+1))
        obs["enrichment"]=obs.observed_w/obs.mean_exp_w.replace(0,np.nan)
        # FDR (Benjamini-Hochberg)
        obs=obs.sort_values("p_emp").reset_index(drop=True)
        m=len(obs)
        obs["fdr"]=obs.p_emp* m/(obs.index+1)
        obs.to_parquet(out_dir/"convergent_kegg.parquet",index=False)
        print(f"✔ convergent_kegg.parquet ({obs.shape[0]} KOs)")

    # ── finalise weighted means & dN/dS ─────────────────────────────────────
    def finalise_means(df,sum_cols):
        for col in sum_cols:
            df[col.replace("sum_","mean_")]=df[col]/df["n_snps"]
            df.drop(columns=[col],inplace=True)
        return df
    for name,df in runs.items():
        if name in ["sweep_core"] or df is None: continue
        if name in ["gene_traj","go_stats","kegg_stats"]:
            df=finalise_means(df,["sum_abs_slope","sum_p","sum_delta"])
        if name.endswith("dnds"):
            df["dnds"]=df.apply(lambda r:r.dn/r.ds if r.ds else np.nan,axis=1)
        pq.write_table(pa.Table.from_pandas(df,preserve_index=False),
                       out_dir/f"{name}.parquet",compression="snappy")
        print(f"✔ {name}.parquet ({df.shape[0]:,} rows)")

    # ── burden per patient/week ────────────────────────────────────────────
    long_ds=ds.dataset(out_dir/"mutation_long",format="parquet")
    burden=(long_ds.to_table(columns=["patient_id","week_num"])
                   .to_pandas()
                   .groupby(["patient_id","week_num"],as_index=False)
                   .size().rename(columns={"size":"n_snvs"}))
    burden.to_parquet(out_dir/"patient_week.parquet",index=False)
    print(f"✔ patient_week.parquet ({burden.shape[0]:,} rows)")
    print("🏁 Finished — all summaries written to", out_dir)

if __name__=="__main__":
    main()