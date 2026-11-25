#!/usr/bin/env bash
# =============================================================================
# server_build_candidate_trees.sh — build candidate ASV trees on server
#
# Usage (server):
#   bash server_build_candidate_trees.sh <workdir> [threads]
#
# Expects: <workdir>/inputs/{region}_{filter}.fa
# Writes:  <workdir>/outputs/{region}_{filter}/{aln.fa,trimmed.fa,tree.nwk,*.log,DONE}
#
# 9 jobs = 3 regions x 3 filters. Runs the three regions in parallel per filter
# (64 threads each by default; server has 256 cores).
# =============================================================================
set -euo pipefail

WORKDIR="${1:?usage: bash server_build_candidate_trees.sh <workdir> [threads]}"
THREADS="${2:-64}"

MAFFT="${MAFFT:-/opt/biosoft/mafft/bin/mafft}"
TRIMAL="${TRIMAL:-${TRIMAL:-trimal}}"
FASTTREE="${FASTTREE:-/opt/biosoft/FastTree/FastTreeMP}"
TRIMAL_MODE="${TRIMAL_MODE:--automated1}"

for f in "$MAFFT" "$TRIMAL" "$FASTTREE"; do
  [ -x "$f" ] || { echo "missing executable: $f" >&2; exit 1; }
done
mkdir -p "$WORKDIR/inputs" "$WORKDIR/outputs"

build_one() {
  local reg="$1" fl="$2"
  local fa="$WORKDIR/inputs/${reg}_${fl}.fa"
  local out="$WORKDIR/outputs/${reg}_${fl}"
  local n
  n=$(grep -c '^>' "$fa")
  mkdir -p "$out"
  echo "[$reg/$fl] input=$fa n=$n start=$(date '+%F %T')"

  if [ ! -s "$out/aln.fa" ]; then
    echo "[$reg/$fl] mafft $(date '+%F %T')"
    "$MAFFT" --auto --thread "$THREADS" "$fa" > "$out/aln.fa" 2> "$out/mafft.log"
  fi
  [ -s "$out/aln.fa" ] || { echo "[$reg/$fl] mafft failed" >&2; exit 1; }

  if [ ! -s "$out/trimmed.fa" ]; then
    echo "[$reg/$fl] trimal $TRIMAL_MODE $(date '+%F %T')"
    "$TRIMAL" -in "$out/aln.fa" -out "$out/trimmed.fa" "$TRIMAL_MODE" 2> "$out/trimal.log"
  fi
  [ -s "$out/trimmed.fa" ] || { echo "[$reg/$fl] trimal failed" >&2; exit 1; }

  if [ ! -s "$out/tree.nwk" ]; then
    echo "[$reg/$fl] FastTreeMP $(date '+%F %T')"
    OMP_NUM_THREADS="$THREADS" "$FASTTREE" -nt -gtr -gamma "$out/trimmed.fa" \
      > "$out/tree.nwk" 2> "$out/fasttree.log"
  fi
  [ -s "$out/tree.nwk" ] || { echo "[$reg/$fl] FastTree failed" >&2; exit 1; }

  echo "$n" > "$out/n_input.txt"
  echo "[$reg/$fl] done=$(date '+%F %T') n_input=$n tree_bytes=$(stat -c%s "$out/tree.nwk")"
  touch "$out/DONE"
}

for fl in main20_5_3 loose10_3_2 strict30_5_3; do
  echo "=== filter $fl $(date '+%F %T') ==="
  for reg in v4 v5_v7 its1; do
    build_one "$reg" "$fl" &
  done
  wait
done

echo "ALL DONE $(date '+%F %T')"
