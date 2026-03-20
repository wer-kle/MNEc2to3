#!/usr/bin/env python3

from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Audit and classify GGPPlus EquCab2->EquCab3 remap results.")
    parser.add_argument(
        "--input",
        default="/Users/weronikaklecel/blast_ggpplus_run/GGPPlus.bed_plus_blast.tsv",
        help="Master merged BED+BLAST table from process_blast_illumina.py"
    )
    parser.add_argument(
        "--outdir",
        default="/Users/weronikaklecel/blast_ggpplus_run/audit",
        help="Output directory"
    )
    parser.add_argument(
        "--prefix",
        default="GGPPlus",
        help="Output file prefix"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    outdir = Path(args.outdir)
    prefix = args.prefix
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", low_memory=False)

    # -----------------------------
    # Defensive type cleanup
    # -----------------------------
    bool_cols = [
        "EquCab2_unique",
        "EquCab3_unique",
        "unique_in_both",
        "bed_chr_matches_EquCab2_blast",
        "EquCab2_exact_match_bed",
        "EquCab2_within1bp_bed",
    ]
    for col in bool_cols:
        if col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].astype(str).str.strip().str.lower().map({
                    "true": True,
                    "false": False,
                    "1": True,
                    "0": False
                })
            df[col] = df[col].fillna(False)

    numeric_cols = [
        "EquCab2_bed_pos",
        "EquCab2_snp_pos",
        "EquCab3_snp_pos",
        "EquCab2_hit_count",
        "EquCab3_hit_count",
        "EquCab2_blast_minus_bed_pos",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    required = [
        "probe_id",
        "EquCab2_unique",
        "EquCab3_unique",
        "unique_in_both",
        "bed_chr_matches_EquCab2_blast",
        "EquCab2_exact_match_bed",
        "EquCab2_within1bp_bed",
        "EquCab2_blast_minus_bed_pos",
        "EquCab2_hit_count",
        "EquCab3_hit_count",
        "EquCab2_bed_chr",
        "EquCab2_bed_pos",
        "EquCab3_chrom",
        "EquCab3_snp_pos",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}")

    # -----------------------------
    # Core categories
    # -----------------------------
    high_conf_exact = df[
        (df["unique_in_both"]) &
        (df["EquCab2_exact_match_bed"])
    ].copy()

    high_conf_off1 = df[
        (df["unique_in_both"]) &
        (~df["EquCab2_exact_match_bed"]) &
        (df["bed_chr_matches_EquCab2_blast"]) &
        (df["EquCab2_blast_minus_bed_pos"].abs() == 1)
    ].copy()

    discordant_unique_both = df[
        (df["unique_in_both"]) &
        (
            (~df["bed_chr_matches_EquCab2_blast"]) |
            (df["EquCab2_blast_minus_bed_pos"].abs() > 1) |
            (df["EquCab2_blast_minus_bed_pos"].isna())
        )
    ].copy()

    unique_ec3_only = df[
        (df["EquCab3_unique"]) &
        (~df["unique_in_both"])
    ].copy()

    unresolved_nonunique_ec3 = df[
        (~df["EquCab3_unique"])
    ].copy()

    # -----------------------------
    # Useful subcategories
    # -----------------------------
    unique_ec3_only_eq2_multihit = unique_ec3_only[
        unique_ec3_only["EquCab2_hit_count"] > 1
    ].copy()

    unique_ec3_only_eq2_nohit = unique_ec3_only[
        unique_ec3_only["EquCab2_hit_count"] == 0
    ].copy()

    unique_ec3_only_eq2_single_but_flag_issue = unique_ec3_only[
        unique_ec3_only["EquCab2_hit_count"] == 1
    ].copy()

    unresolved_ec3_multihit = unresolved_nonunique_ec3[
        unresolved_nonunique_ec3["EquCab3_hit_count"] > 1
    ].copy()

    unresolved_ec3_nohit = unresolved_nonunique_ec3[
        unresolved_nonunique_ec3["EquCab3_hit_count"] == 0
    ].copy()

    # -----------------------------
    # Compact deliverable tables
    # -----------------------------
    deliverable_cols = [
        "probe_id",
        "EquCab2_bed_chr",
        "EquCab2_bed_pos",
        "EquCab2_hit_count",
        "EquCab3_hit_count",
        "EquCab2_unique",
        "EquCab3_unique",
        "unique_in_both",
        "EquCab2_chrom",
        "EquCab2_snp_pos",
        "EquCab2_strand",
        "EquCab3_chrom",
        "EquCab3_snp_pos",
        "EquCab3_strand",
        "bed_chr_matches_EquCab2_blast",
        "EquCab2_blast_minus_bed_pos",
        "EquCab2_exact_match_bed",
        "EquCab2_within1bp_bed",
    ]
    deliverable_cols = [c for c in deliverable_cols if c in df.columns]

    high_conf_exact = high_conf_exact[deliverable_cols].sort_values("probe_id")
    high_conf_off1 = high_conf_off1[deliverable_cols].sort_values("probe_id")
    discordant_unique_both = discordant_unique_both[deliverable_cols].sort_values("probe_id")
    unique_ec3_only = unique_ec3_only[deliverable_cols].sort_values("probe_id")
    unresolved_nonunique_ec3 = unresolved_nonunique_ec3[deliverable_cols].sort_values("probe_id")

    unique_ec3_only_eq2_multihit = unique_ec3_only_eq2_multihit[deliverable_cols].sort_values("probe_id")
    unique_ec3_only_eq2_nohit = unique_ec3_only_eq2_nohit[deliverable_cols].sort_values("probe_id")
    unique_ec3_only_eq2_single_but_flag_issue = unique_ec3_only_eq2_single_but_flag_issue[deliverable_cols].sort_values("probe_id")
    unresolved_ec3_multihit = unresolved_ec3_multihit[deliverable_cols].sort_values("probe_id")
    unresolved_ec3_nohit = unresolved_ec3_nohit[deliverable_cols].sort_values("probe_id")

    # -----------------------------
    # Write tables
    # -----------------------------
    files = {
        f"{prefix}.audit.high_confidence_exact.tsv": high_conf_exact,
        f"{prefix}.audit.high_confidence_off_by_1.tsv": high_conf_off1,
        f"{prefix}.audit.discordant_unique_both.tsv": discordant_unique_both,
        f"{prefix}.audit.unique_EquCab3_only.tsv": unique_ec3_only,
        f"{prefix}.audit.unresolved_nonunique_EquCab3.tsv": unresolved_nonunique_ec3,
        f"{prefix}.audit.unique_EquCab3_only.eq2_multihit.tsv": unique_ec3_only_eq2_multihit,
        f"{prefix}.audit.unique_EquCab3_only.eq2_nohit.tsv": unique_ec3_only_eq2_nohit,
        f"{prefix}.audit.unique_EquCab3_only.eq2_single_but_flag_issue.tsv": unique_ec3_only_eq2_single_but_flag_issue,
        f"{prefix}.audit.unresolved_EquCab3_multihit.tsv": unresolved_ec3_multihit,
        f"{prefix}.audit.unresolved_EquCab3_nohit.tsv": unresolved_ec3_nohit,
    }

    for fname, table in files.items():
        table.to_csv(outdir / fname, sep="\t", index=False)

    # -----------------------------
    # Summary
    # -----------------------------
    n_total = len(df)
    n_exact = len(high_conf_exact)
    n_off1 = len(high_conf_off1)
    n_discord = len(discordant_unique_both)
    n_ec3_only = len(unique_ec3_only)
    n_unresolved = len(unresolved_nonunique_ec3)

    with open(outdir / f"{prefix}.audit.summary.txt", "w") as out:
        out.write("GGPPlus remap audit summary\n")
        out.write("===========================\n\n")

        out.write(f"Total probes in master table: {n_total}\n\n")

        out.write("Primary classes\n")
        out.write("---------------\n")
        out.write(f"High-confidence exact:              {n_exact}\n")
        out.write(f"High-confidence off-by-1:           {n_off1}\n")
        out.write(f"Discordant but unique in both:      {n_discord}\n")
        out.write(f"Unique in EquCab3 only:             {n_ec3_only}\n")
        out.write(f"Unresolved / non-unique in EquCab3: {n_unresolved}\n\n")

        out.write("Consistency checks\n")
        out.write("------------------\n")
        out.write(f"Exact + off-by-1 + discordant unique-both: {n_exact + n_off1 + n_discord}\n")
        out.write(f"Total unique in both expected from upstream: {(df['unique_in_both']).sum()}\n\n")

        out.write("Useful subcategories\n")
        out.write("--------------------\n")
        out.write(f"Unique in EquCab3 only, EquCab2 multihit: {len(unique_ec3_only_eq2_multihit)}\n")
        out.write(f"Unique in EquCab3 only, EquCab2 no hit:   {len(unique_ec3_only_eq2_nohit)}\n")
        out.write(f"Unique in EquCab3 only, EquCab2 single:   {len(unique_ec3_only_eq2_single_but_flag_issue)}\n")
        out.write(f"Unresolved EquCab3 multihit:              {len(unresolved_ec3_multihit)}\n")
        out.write(f"Unresolved EquCab3 no hit:                {len(unresolved_ec3_nohit)}\n\n")

        out.write("Suggested use\n")
        out.write("-------------\n")
        out.write("Primary high-confidence remap set: use high_confidence_exact.tsv\n")
        out.write("Secondary acceptable set after review: high_confidence_off_by_1.tsv\n")
        out.write("Manual review set: discordant_unique_both.tsv\n")
        out.write("Do not use for strict remap without extra justification: unique_EquCab3_only.tsv and unresolved_nonunique_EquCab3.tsv\n\n")

        out.write("Files written\n")
        out.write("-------------\n")
        for fname in files:
            out.write(f"{outdir / fname}\n")
        out.write(f"{outdir / f'{prefix}.audit.summary.txt'}\n")

    print("Done.")
    print(f"Wrote audit files to: {outdir}")
    print(f"Summary: {outdir / f'{prefix}.audit.summary.txt'}")


if __name__ == "__main__":
    main()