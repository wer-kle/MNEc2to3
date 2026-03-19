#!/usr/bin/env bash
set -euo pipefail

QUERY_FASTA="/Volumes/Public/SNP_array_data/horse_data/ref/manifest_files/GGPPlus_manifest_AlleleA_ProbeSeq.fasta"
REF_EQCAB3="/Volumes/Public/SNP_array_data/horse_data/ref/reference_genomes/GCF_002863925.1_EquCab3.0_genomic.fna"
REF_EQCAB2="/Volumes/Public/SNP_array_data/horse_data/ref/reference_genomes/GCF_000002305.2_EquCab2.0_genomic.fna"

WORKDIR="$HOME/blast_ggpplus_run"
THREADS=4

mkdir -p "$WORKDIR"
cd "$WORKDIR"

echo "Building local BLAST databases (version 4)..."
makeblastdb -in "$REF_EQCAB2" -dbtype nucl -parse_seqids -blastdb_version 4 -out EquCab2
makeblastdb -in "$REF_EQCAB3" -dbtype nucl -parse_seqids -blastdb_version 4 -out EquCab3

echo "Running BLAST..."
for DB in EquCab2 EquCab3; do
  blastn \
    -task blastn \
    -db "$DB" \
    -query "$QUERY_FASTA" \
    -word_size 50 \
    -perc_identity 100 \
    -dust no \
    -soft_masking false \
    -strand both \
    -num_threads "$THREADS" \
    -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send sstrand" \
    -out "GGPPlus.raw.${DB}.tsv"

  awk 'BEGIN{OFS="\t"} $3==100 && $4==$5 && $6==1 && $7==$5' \
    "GGPPlus.raw.${DB}.tsv" > "GGPPlus.exact.${DB}.tsv"
done

echo "Done."
echo "Results are in: $WORKDIR"
