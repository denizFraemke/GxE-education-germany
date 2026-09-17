#!/usr/bin/env bash
# =============================================================================
# format_height_sumstats.sh — Build per-chromosome SBayesR .ma inputs for Height
#                             from the raw GIANT/Yengo 2022 GWAS (a share file)
# =============================================================================
# Produces the 22 files the pipeline expects:
#   ${SUMSTATS_DIR}/Height_GIANT/HEIGHT_SBayesRformat.chr{1..22}.txt
#
# WHY THIS EXISTS
# ---------------
# The pipeline is reproducible from the data-share inputs: every other trait
# ships a correctly pre-formatted per-chromosome set on the share
# (013-PGS/data/formatted/<TRAIT>/). Height is the exception: that pre-formatted
# copy is MALFORMED — it keeps GIANT's leading CHR column (9 columns instead of
# 8) and does not drop SNPs with missing effect estimates (NA BETA/SE/P). With
# GCTB >= 2.5 that shifts every column by one → "0 matched SNPs" and, on the
# chromosome that hits an NA row first, a segfault. Height is therefore built
# from the RAW GIANT file, which is clean and on the share.
#
# RAW INPUT (GIANT/Yengo 2022, 11 tab-separated columns):
#   SNPID  RSID  CHR  POS  EFFECT_ALLELE  OTHER_ALLELE  EFFECT_ALLELE_FREQ
#   BETA  SE  P  N
#   data share: Projects/03_data/009_SUMSTATS/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz
#
# OUTPUT (standard SBayesR/GCTA .ma, 8 columns, header + data, one per chr):
#   SNP  A1  A2  freq  b  se  p  N
#   - SNP = RSID  (the UKB-50k LD reference is rsID-keyed; the other traits'
#                  .ma files are rsID-keyed too, which is why they match)
#   - CHR is used ONLY to split into per-chromosome files; it is NOT written
#   - autosomes 1-22 only; rows missing any of freq/BETA/SE/P/N are dropped
#
# USAGE
#   bash scripts/format_height_sumstats.sh [RAW_GIANT_FILE]
#   # RAW_GIANT_FILE defaults to $GIANT_HEIGHT_RAW (set it in user_config.sh).
#   # Run it on any machine that can read the raw GIANT file (the Mac with the
#   # share mounted, or the cluster after staging the raw file), then the resulting
#   # Height_GIANT/ folder is staged/used like any other trait.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config.sh"

RAW="${1:-${GIANT_HEIGHT_RAW:-}}"
OUT_DIR="${SUMSTATS_DIR}/Height_GIANT"

if [[ -z "$RAW" || ! -f "$RAW" ]]; then
    cat >&2 <<EOF
ERROR: raw GIANT height file not found: '${RAW:-<unset>}'
Pass it as the first argument, or set GIANT_HEIGHT_RAW in user_config.sh.
On the data share it is:
  Projects/03_data/009_SUMSTATS/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz
EOF
    exit 1
fi

mkdir -p "$OUT_DIR"
read_raw() { case "$RAW" in *.gz) gzip -dc -- "$RAW";; *) cat -- "$RAW";; esac; }

TMP="$(mktemp)"; trap 'rm -f "$TMP"' EXIT

echo "Formatting Height (raw -> per-chr SBayesR .ma)"
echo "  raw : $RAW"
echo "  out : $OUT_DIR"

# Pass 1: locate columns BY NAME (robust to column re-ordering), keep autosomes,
# drop rows with any missing stat, emit "chr \t SNP A1 A2 freq b se p N".
read_raw | awk -F'\t' '
NR==1 {
  for (i=1; i<=NF; i++) col[$i] = i
  split("RSID CHR EFFECT_ALLELE OTHER_ALLELE EFFECT_ALLELE_FREQ BETA SE P N", need, " ")
  for (k in need) if (!(need[k] in col)) {
    print "ERROR: expected GIANT column not found: " need[k] > "/dev/stderr"; exit 2
  }
  next
}
{
  c = $(col["CHR"])
  if (c !~ /^[0-9]+$/ || c+0 < 1 || c+0 > 22) next
  fr=$(col["EFFECT_ALLELE_FREQ"]); b=$(col["BETA"]); se=$(col["SE"]); p=$(col["P"]); n=$(col["N"])
  if (fr==""||fr=="NA"||b==""||b=="NA"||se==""||se=="NA"||p==""||p=="NA"||n==""||n=="NA") { drop++; next }
  print c "\t" $(col["RSID"]) "\t" $(col["EFFECT_ALLELE"]) "\t" $(col["OTHER_ALLELE"]) "\t" fr "\t" b "\t" se "\t" p "\t" n
  keep++
}
END { printf "  kept %d SNPs, dropped %d (missing stats)\n", keep+0, drop+0 > "/dev/stderr" }
' > "$TMP"

# Pass 2: one 8-column file per chromosome (strip the leading chr field), header first.
for c in $(seq 1 22); do
  out="${OUT_DIR}/HEIGHT_SBayesRformat.chr${c}.txt"
  { printf 'SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n'
    awk -F'\t' -v c="$c" '$1==c { sub(/^[0-9]+\t/, ""); print }' "$TMP"
  } > "$out"
done

# Verify: every file has exactly 8 columns.
for c in $(seq 1 22); do
  f="${OUT_DIR}/HEIGHT_SBayesRformat.chr${c}.txt"
  nf=$(awk -F'\t' 'NR==1{print NF; exit}' "$f")
  [[ "$nf" == "8" ]] || { echo "ERROR: ${f} has ${nf} columns (expected 8)" >&2; exit 3; }
done

echo "  OK: 22 per-chromosome files written, 8 columns each."
echo "  header: $(head -1 "${OUT_DIR}/HEIGHT_SBayesRformat.chr1.txt")"
