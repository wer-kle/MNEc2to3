#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


BLAST_COLS = [
    "probe_id",
    "sseqid",
    "pident",
    "align_len",
    "qlen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstrand",
]


def normalize_chr_label(x):
    if pd.isna(x):
        return pd.NA

    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return pd.NA

    s = s.replace("chromosome ", "").replace("Chromosome ", "")
    s = re.sub(r"^chr", "", s, flags=re.IGNORECASE).strip()

    if re.fullmatch(r"\d+", s):
        return f"chr{int(s)}"
    if s.upper() in {"X", "Y"}:
        return f"chr{s.upper()}"
    if s.upper() in {"M", "MT", "MITO"}:
        return "chrM"

    return str(x).strip()


def parse_accession_from_sseqid(sseqid: str):
    if pd.isna(sseqid):
        return None

    s = str(sseqid)
    m = re.search(r"([A-Z]{2}_[0-9]+\.[0-9]+)", s)
    if m:
        return m.group(1)
    return None


def build_accession_map_from_fasta(fasta_path: Path):
    mapping = {}

    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            seqid = header.split()[0]

            accession = parse_accession_from_sseqid(seqid)
            if accession is None:
                accession = parse_accession_from_sseqid(header)

            chrom_label = None

            m_chr = re.search(r"chromosome\s+([0-9XY]+)\b", header, flags=re.IGNORECASE)
            if m_chr:
                chrom = m_chr.group(1).upper()
                chrom_label = f"chr{chrom}"
            elif re.search(r"mitochond", header, flags=re.IGNORECASE):
                chrom_label = "chrM"

            if accession and chrom_label:
                mapping[accession] = chrom_label
                mapping[f"ref|{accession}|"] = chrom_label

            if seqid and chrom_label:
                mapping[seqid] = chrom_label

    return mapping


def normalize_blast_chr(sseqid, accession_map):
    if pd.isna(sseqid):
        return pd.NA

    s = str(sseqid).strip()

    if s in accession_map:
        return accession_map[s]

    acc = parse_accession_from_sseqid(s)
    if acc and acc in accession_map:
        return accession_map[acc]

    return normalize_chr_label(s)


def load_bed4(bed_path: Path):
    bed = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["bed_chr_raw", "bed_start_0based", "bed_end_1based", "probe_id"],
        dtype={0: str, 1: int, 2: int, 3: str},
    )

    bed["probe_id"] = bed["probe_id"].astype(str).str.strip()
    bed["EquCab2_bed_chr"] = bed["bed_chr_raw"].apply(normalize_chr_label)

    # BED is 0-based start, half-open end. For SNP BED4, 1-based SNP position is start + 1.
    bed["EquCab2_bed_pos"] = bed["bed_start_0based"] + 1

    return bed


def load_exact_blast(path: Path, accession_map, genome_name: str):
    df = pd.read_csv(path, sep="\t", header=None, names=BLAST_COLS)

    df["probe_id"] = df["probe_id"].astype(str).str.strip()
    df["chrom_raw"] = df["sseqid"].astype(str)
    df["chrom"] = df["sseqid"].apply(lambda x: normalize_blast_chr(x, accession_map))

    df["probe_start_1based"] = df[["sstart", "send"]].min(axis=1)
    df["probe_end_1based"] = df[["sstart", "send"]].max(axis=1)

    # Two competing models:
    # A: SNP lies just beyond probe 3' end
    # B: SNP lies just before probe 5' end
    df["snp_pos_model_A"] = np.where(
        df["sstrand"].eq("plus"),
        df["probe_end_1based"] + 1,
        df["probe_start_1based"] - 1,
    )

    df["snp_pos_model_B"] = np.where(
        df["sstrand"].eq("plus"),
        df["probe_start_1based"] - 1,
        df["probe_end_1based"] + 1,
    )

    counts = df.groupby("probe_id").size().reset_index(name=f"{genome_name}_hit_count")
    df = df.merge(counts, on="probe_id", how="left")
    df[f"{genome_name}_unique"] = df[f"{genome_name}_hit_count"] == 1

    return df


def summarize_model_fit(eq2_unique_with_bed: pd.DataFrame, pos_col: str):
    same_chr = eq2_unique_with_bed["EquCab2_bed_chr"] == eq2_unique_with_bed["chrom"]
    pos_diff = eq2_unique_with_bed[pos_col] - eq2_unique_with_bed["EquCab2_bed_pos"]

    exact = (same_chr & (pos_diff == 0)).sum()
    within1 = (same_chr & (pos_diff.abs() <= 1)).sum()
    within5 = (same_chr & (pos_diff.abs() <= 5)).sum()

    median_abs_diff = pos_diff[same_chr].abs().median() if same_chr.any() else np.nan

    return {
        "exact": int(exact),
        "within1": int(within1),
        "within5": int(within5),
        "median_abs_diff_same_chr": float(median_abs_diff) if pd.notna(median_abs_diff) else np.nan,
    }


def choose_best_model(eq2_unique_with_bed: pd.DataFrame):
    score_a = summarize_model_fit(eq2_unique_with_bed, "snp_pos_model_A")
    score_b = summarize_model_fit(eq2_unique_with_bed, "snp_pos_model_B")

    a_key = (
        score_a["exact"],
        score_a["within1"],
        -score_a["median_abs_diff_same_chr"] if pd.notna(score_a["median_abs_diff_same_chr"]) else -1e15,
    )
    b_key = (
        score_b["exact"],
        score_b["within1"],
        -score_b["median_abs_diff_same_chr"] if pd.notna(score_b["median_abs_diff_same_chr"]) else -1e15,
    )

    chosen = "A" if a_key >= b_key else "B"
    return chosen, score_a, score_b


def make_unique_table(df: pd.DataFrame, genome_name: str, chosen_model: str):
    pos_col = f"snp_pos_model_{chosen_model}"
    unique_df = df[df[f"{genome_name}_unique"]].copy()

    out = unique_df[
        [
            "probe_id",
            "chrom_raw",
            "chrom",
            "probe_start_1based",
            "probe_end_1based",
            "sstart",
            "send",
            "sstrand",
            pos_col,
        ]
    ].copy()

    out = out.rename(columns={
        "chrom_raw": f"{genome_name}_chrom_raw",
        "chrom": f"{genome_name}_chrom",
        "probe_start_1based": f"{genome_name}_probe_start_1based",
        "probe_end_1based": f"{genome_name}_probe_end_1based",
        "sstart": f"{genome_name}_sstart_raw",
        "send": f"{genome_name}_send_raw",
        "sstrand": f"{genome_name}_strand",
        pos_col: f"{genome_name}_snp_pos",
    })

    return out.sort_values("probe_id").reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser(description="GGPPlus remap using BED4 EquCab2 coordinates and BLAST exact hits.")
    parser.add_argument(
        "--bed4",
        default="/Volumes/wklecel/milestone1/Univ_of_Calgary_Poissant_EQNG75V03_20181018_SNP_Map.bed",
        help="BED4 file with authoritative EquCab2 coordinates."
    )
    parser.add_argument(
        "--blast_dir",
        default="/Users/weronikaklecel/blast_ggpplus_run",
        help="Directory containing GGPPlus.exact.EquCab2.tsv and GGPPlus.exact.EquCab3.tsv"
    )
    parser.add_argument(
        "--ref_eqcab2",
        default="/Volumes/Public/SNP_array_data/horse_data/ref/reference_genomes/GCF_000002305.2_EquCab2.0_genomic.fna",
        help="EquCab2 reference FASTA"
    )
    parser.add_argument(
        "--ref_eqcab3",
        default="/Volumes/Public/SNP_array_data/horse_data/ref/reference_genomes/GCF_002863925.1_EquCab3.0_genomic.fna",
        help="EquCab3 reference FASTA"
    )
    parser.add_argument(
        "--prefix",
        default="GGPPlus",
        help="Prefix for output files"
    )
    args = parser.parse_args()

    bed_path = Path(args.bed4)
    blast_dir = Path(args.blast_dir)
    ref_eqcab2 = Path(args.ref_eqcab2)
    ref_eqcab3 = Path(args.ref_eqcab3)
    prefix = args.prefix

    eq2_exact = blast_dir / f"{prefix}.exact.EquCab2.tsv"
    eq3_exact = blast_dir / f"{prefix}.exact.EquCab3.tsv"

    out_counts = blast_dir / f"{prefix}.probe_blast_counts.tsv"
    out_eq2_exact_norm = blast_dir / f"{prefix}.exact.EquCab2.normalized.tsv"
    out_eq3_exact_norm = blast_dir / f"{prefix}.exact.EquCab3.normalized.tsv"
    out_eq2_unique = blast_dir / f"{prefix}.unique_hits.EquCab2.normalized.tsv"
    out_eq3_unique = blast_dir / f"{prefix}.unique_hits.EquCab3.normalized.tsv"
    out_merged = blast_dir / f"{prefix}.bed_plus_blast.tsv"
    out_remap_ec3 = blast_dir / f"{prefix}.remap_to_EquCab3.unique_ec3.tsv"
    out_remap_both = blast_dir / f"{prefix}.remap_to_EquCab3.unique_both.tsv"
    out_summary = blast_dir / f"{prefix}.remap_summary.txt"

    bed = load_bed4(bed_path)

    eq2_map = build_accession_map_from_fasta(ref_eqcab2)
    eq3_map = build_accession_map_from_fasta(ref_eqcab3)

    eq2 = load_exact_blast(eq2_exact, eq2_map, "EquCab2")
    eq3 = load_exact_blast(eq3_exact, eq3_map, "EquCab3")

    eq2[
        [
            "probe_id", "chrom_raw", "chrom", "pident", "align_len", "qlen",
            "qstart", "qend", "sstart", "send", "sstrand",
            "probe_start_1based", "probe_end_1based",
            "snp_pos_model_A", "snp_pos_model_B", "EquCab2_hit_count", "EquCab2_unique"
        ]
    ].to_csv(out_eq2_exact_norm, sep="\t", index=False)

    eq3[
        [
            "probe_id", "chrom_raw", "chrom", "pident", "align_len", "qlen",
            "qstart", "qend", "sstart", "send", "sstrand",
            "probe_start_1based", "probe_end_1based",
            "snp_pos_model_A", "snp_pos_model_B", "EquCab3_hit_count", "EquCab3_unique"
        ]
    ].to_csv(out_eq3_exact_norm, sep="\t", index=False)

    eq2_unique = eq2[eq2["EquCab2_unique"]].copy()
    eq2_unique_bed = eq2_unique.merge(
        bed[["probe_id", "EquCab2_bed_chr", "EquCab2_bed_pos"]],
        on="probe_id",
        how="inner"
    )

    chosen_model, score_a, score_b = choose_best_model(eq2_unique_bed)

    eq2_unique_out = make_unique_table(eq2, "EquCab2", chosen_model)
    eq3_unique_out = make_unique_table(eq3, "EquCab3", chosen_model)

    eq2_unique_out.to_csv(out_eq2_unique, sep="\t", index=False)
    eq3_unique_out.to_csv(out_eq3_unique, sep="\t", index=False)

    all_probe_ids = sorted(set(eq2["probe_id"]).union(set(eq3["probe_id"])).union(set(bed["probe_id"])))

    eq2_counts = eq2.groupby("probe_id").size().reset_index(name="EquCab2_hit_count")
    eq3_counts = eq3.groupby("probe_id").size().reset_index(name="EquCab3_hit_count")

    counts = pd.DataFrame({"probe_id": all_probe_ids})
    counts = counts.merge(eq2_counts, on="probe_id", how="left")
    counts = counts.merge(eq3_counts, on="probe_id", how="left")
    counts["EquCab2_hit_count"] = counts["EquCab2_hit_count"].fillna(0).astype(int)
    counts["EquCab3_hit_count"] = counts["EquCab3_hit_count"].fillna(0).astype(int)
    counts["EquCab2_unique"] = counts["EquCab2_hit_count"] == 1
    counts["EquCab3_unique"] = counts["EquCab3_hit_count"] == 1
    counts["unique_in_both"] = counts["EquCab2_unique"] & counts["EquCab3_unique"]

    counts.to_csv(out_counts, sep="\t", index=False)

    merged = bed.merge(counts, on="probe_id", how="left")
    merged = merged.merge(eq2_unique_out, on="probe_id", how="left")
    merged = merged.merge(eq3_unique_out, on="probe_id", how="left")

    required_cols = [
        "EquCab2_hit_count", "EquCab3_hit_count",
        "EquCab2_unique", "EquCab3_unique", "unique_in_both"
    ]
    missing = [c for c in required_cols if c not in merged.columns]
    if missing:
        raise RuntimeError(f"Missing expected columns after merge: {missing}")

    merged["bed_chr_matches_EquCab2_blast"] = merged["EquCab2_bed_chr"] == merged["EquCab2_chrom"]
    merged["EquCab2_blast_minus_bed_pos"] = merged["EquCab2_snp_pos"] - merged["EquCab2_bed_pos"]
    merged["EquCab2_exact_match_bed"] = (
        merged["bed_chr_matches_EquCab2_blast"] &
        (merged["EquCab2_blast_minus_bed_pos"] == 0)
    )
    merged["EquCab2_within1bp_bed"] = (
        merged["bed_chr_matches_EquCab2_blast"] &
        (merged["EquCab2_blast_minus_bed_pos"].abs() <= 1)
    )

    merged.to_csv(out_merged, sep="\t", index=False)

    keep_cols = [
        "probe_id",
        "bed_chr_raw",
        "EquCab2_bed_chr",
        "bed_start_0based",
        "bed_end_1based",
        "EquCab2_bed_pos",
        "EquCab2_hit_count",
        "EquCab3_hit_count",
        "EquCab2_unique",
        "EquCab3_unique",
        "unique_in_both",
        "EquCab2_chrom_raw",
        "EquCab2_chrom",
        "EquCab2_snp_pos",
        "EquCab2_strand",
        "EquCab3_chrom_raw",
        "EquCab3_chrom",
        "EquCab3_snp_pos",
        "EquCab3_strand",
        "bed_chr_matches_EquCab2_blast",
        "EquCab2_blast_minus_bed_pos",
        "EquCab2_exact_match_bed",
        "EquCab2_within1bp_bed",
    ]

    remap_ec3 = merged[merged["EquCab3_unique"]].copy()
    remap_both = merged[merged["unique_in_both"]].copy()

    remap_ec3[keep_cols].to_csv(out_remap_ec3, sep="\t", index=False)
    remap_both[keep_cols].to_csv(out_remap_both, sep="\t", index=False)

    with open(out_summary, "w") as out:
        out.write("GGPPlus BED-aware remap summary\n")
        out.write("==============================\n\n")
        out.write(f"Chosen SNP-position model: {chosen_model}\n")
        out.write("Model A = SNP adjacent to probe 3' end\n")
        out.write("Model B = SNP adjacent to probe 5' end\n\n")

        out.write("Validation against EquCab2 BED coordinates:\n")
        out.write(f"  Model A exact matches: {score_a['exact']}\n")
        out.write(f"  Model A within 1 bp:   {score_a['within1']}\n")
        out.write(f"  Model A within 5 bp:   {score_a['within5']}\n")
        out.write(f"  Model A median abs diff on same chr: {score_a['median_abs_diff_same_chr']}\n\n")

        out.write(f"  Model B exact matches: {score_b['exact']}\n")
        out.write(f"  Model B within 1 bp:   {score_b['within1']}\n")
        out.write(f"  Model B within 5 bp:   {score_b['within5']}\n")
        out.write(f"  Model B median abs diff on same chr: {score_b['median_abs_diff_same_chr']}\n\n")

        out.write(f"BED rows: {len(bed)}\n")
        out.write(f"EquCab2 exact-hit alignments: {len(eq2)}\n")
        out.write(f"EquCab3 exact-hit alignments: {len(eq3)}\n")
        out.write(f"Probes unique in EquCab2: {(counts['EquCab2_unique']).sum()}\n")
        out.write(f"Probes unique in EquCab3: {(counts['EquCab3_unique']).sum()}\n")
        out.write(f"Probes unique in both:    {(counts['unique_in_both']).sum()}\n\n")

        out.write("Output files:\n")
        out.write(f"  {out_counts}\n")
        out.write(f"  {out_eq2_exact_norm}\n")
        out.write(f"  {out_eq3_exact_norm}\n")
        out.write(f"  {out_eq2_unique}\n")
        out.write(f"  {out_eq3_unique}\n")
        out.write(f"  {out_merged}\n")
        out.write(f"  {out_remap_ec3}\n")
        out.write(f"  {out_remap_both}\n")

    print("Done.")
    print(f"Chosen model: {chosen_model}")
    print(f"Wrote: {out_counts}")
    print(f"Wrote: {out_eq2_exact_norm}")
    print(f"Wrote: {out_eq3_exact_norm}")
    print(f"Wrote: {out_eq2_unique}")
    print(f"Wrote: {out_eq3_unique}")
    print(f"Wrote: {out_merged}")
    print(f"Wrote: {out_remap_ec3}")
    print(f"Wrote: {out_remap_both}")
    print(f"Wrote: {out_summary}")


if __name__ == "__main__":
    main()
