#!/bin/bash

if [ -z "$1" ]; then
	echo "Usage: dxtopology.sh <MAP> <NRES> <BASENAME> <OUTDIR>"
	exit 0
fi

SCRIPTSPATH=$(dirname -- "$( readlink -f -- "$0"; )")
OUTDIR=$4
BASENAME=$3
MAP=$1
NRES=$2
WORKDIR=$(pwd)

cp ${SCRIPTSPATH}/dxtraces.sh ${OUTDIR}/
chmod u+x ${OUTDIR}/dxtraces.sh

#initialize trace
cd $OUTDIR
./dxtraces.sh $MAP $NRES > log
mv traces/refined-1.pdb ${BASENAME}-init-trace.pdb
rm dxtraces.sh
rm -r traces/

cd ${SCRIPTSPATH}
./pdb2dat.sh $OUTDIR/${BASENAME}-init-trace.pdb $OUTDIR/${BASENAME}-init-trace.dat

#estimate weights
python extract-dist-pred.py $OUTDIR/seq.npz $BASENAME $OUTDIR $NRES >& logpy
python calcmat.py $BASENAME $OUTDIR >& logpy 
python counts_parallel.py $BASENAME $OUTDIR >& logpy

#optimize weights
./assign -pred_dist $OUTDIR/dist-pred-${BASENAME}.dat -beads_dist $OUTDIR/dist-calc-${BASENAME}.dat -startweights ${OUTDIR}/counts.txt -outdir ${OUTDIR}/assignment > logassign

#transform weights to trace
python indmax2dat.py $BASENAME $OUTDIR
./coords2pdb_n -i ${OUTDIR}/out-filled-${BASENAME}.dat -o  ${OUTDIR}/out-filled-${BASENAME}.pdb -n ${OUTDIR}/n_filled.txt
./coords2pdb_n -i ${OUTDIR}/out-cleaned-${BASENAME}.dat -o  ${OUTDIR}/out-cleaned-${BASENAME}.pdb -n ${OUTDIR}/n_clean.txt
./coords2pdb_n -i ${OUTDIR}/out-smoothed-${BASENAME}.dat -o  ${OUTDIR}/out-smoothed-${BASENAME}.pdb -n ${OUTDIR}/n_smooth.txt
./coords2pdb_n -i ${OUTDIR}/out-sorted-${BASENAME}.dat -o  ${OUTDIR}/out-sorted-${BASENAME}.pdb -n ${OUTDIR}/n_sorted.txt
./coords2pdb_n -i ${OUTDIR}/out-raw-${BASENAME}.dat -o  ${OUTDIR}/out-raw-${BASENAME}.pdb -n ${OUTDIR}/n_raw.txt

cd $WORKDIR
