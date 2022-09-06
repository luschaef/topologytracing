#!/bin/bash

OUTDIR="ribo-I-4A"
BASENAME="ribo-I-4A"
NRES="125"

#calculate start weights
./pdb2dat.sh $OUTDIR/${BASENAME}-init-trace.pdb $OUTDIR/${BASENAME}-init-trace.dat
#python3 extract-dist-pred.py $OUTDIR/seq.npz $BASENAME $OUTDIR $NRES
python3 calcmat.py $BASENAME $OUTDIR 
python3 counts_parallel.py $BASENAME $OUTDIR

