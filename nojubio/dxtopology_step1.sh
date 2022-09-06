#!/bin/bash

OUTDIR="ribo-I-4A"
BASENAME="ribo-I-4A"
MAP="ribo-I-4A.mrc"
NRES="125"

#mkdir $OUTDIR
cp dxtraces-3.sh ${OUTDIR}/
chmod u+x ${OUTDIR}/dxtraces-3.sh
#initialize trace

cd /local/jubio/luisa/topologytracing/$OUTDIR
./dxtraces-3.sh $MAP $NRES $OUTDIR > log
mv traces/refined-1.pdb ${BASENAME}-init-trace.pdb
rm ${OUTDIR}/dxtraces-3.sh

