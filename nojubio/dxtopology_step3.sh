#!/bin/bash

OUTDIR="26swo35"
BASENAME="26swo35"

#run assignment
cd /local/jubio/luisa/topologytracing/
./assign_clean -pred_dist $OUTDIR/dist-pred-${BASENAME}.dat -beads_dist $OUTDIR/dist-calc-${BASENAME}.dat -startweights ${OUTDIR}/counts.txt -outdir ${OUTDIR}/assignmenttimetestw35tothe3 > logassign


