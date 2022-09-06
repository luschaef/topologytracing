#!/bin/bash

OUTDIR="TMV_2"
BASENAME="tmv"

python3 indmax2dat_new.py $BASENAME $OUTDIR
./coords2pdb -i ${OUTDIR}/out-filled-${BASENAME}.dat -o  ${OUTDIR}/out-filled-${BASENAME}.pdb -n `cat ${OUTDIR}/n_filled.txt`
./coords2pdb -i ${OUTDIR}/out-cleaned-${BASENAME}.dat -o  ${OUTDIR}/out-cleaned-${BASENAME}.pdb -n `cat ${OUTDIR}/n_clean.txt`
./coords2pdb -i ${OUTDIR}/out-smoothed-${BASENAME}.dat -o  ${OUTDIR}/out-smoothed-${BASENAME}.pdb -n `cat ${OUTDIR}/n_smooth.txt`
./coords2pdb -i ${OUTDIR}/out-sorted-${BASENAME}.dat -o  ${OUTDIR}/out-sorted-${BASENAME}.pdb -n `cat ${OUTDIR}/n_sorted.txt`
./coords2pdb -i ${OUTDIR}/out-raw-${BASENAME}.dat -o  ${OUTDIR}/out-raw-${BASENAME}.pdb -n `cat ${OUTDIR}/n_raw.txt`
#rm ${OUTDIR}/n_filled.txt
#rm ${OUTDIR}/n_clean.txt
#rm ${OUTDIR}/n_smooth.txt
#rm ${OUTDIR}/n_sorted.txt

