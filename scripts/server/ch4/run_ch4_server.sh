#!/bin/bash
set -x
cd <SERVER_WORKDIR> || exit 1
export BASE_DIR=$PWD
export GOS_TREE_DIR=$PWD/output
export OPENBLAS_NUM_THREADS=64
export OMP_NUM_THREADS=64
export FAST_LAMBDA=1
RS=${RSCRIPT:-Rscript}
$RS scripts/R/gos_allASV_phylosig_fast.R > output/phylosig_all.log 2>&1
echo "phylosig all exit $?"
SUBSET=positive $RS scripts/R/gos_allASV_phylosig_fast.R > output/phylosig_positive.log 2>&1
echo "phylosig positive exit $?"
GOS_TREE_DIR=$PWD/output TREE_FILE=gos_all_hostn_16s_midpoint.nwk ANN_FILE=annotation_16s_allASV.tsv CAT_COL=category GEN_VALUES=generalist SUBSET=all OUT_FILE=phyloD_generalist_all.tsv D_PERMUT=1000 $RS scripts/R/gos_generalist_D.R > output/D_all.log 2>&1
echo "D all exit $?"
GOS_TREE_DIR=$PWD/output TREE_FILE=gos_all_hostn_16s_midpoint.nwk ANN_FILE=annotation_16s_allASV.tsv CAT_COL=category GEN_VALUES=generalist SUBSET=positive OUT_FILE=phyloD_generalist_positive.tsv D_PERMUT=1000 $RS scripts/R/gos_generalist_D.R > output/D_positive.log 2>&1
echo "D positive exit $?"
echo ALLDONE
