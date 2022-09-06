#!/bin/bash

OUTDIR="26swo35"
BASENAME="26swo35"

#run assignment
cat << EOF > jubio-assign-${BASENAME}.sh
#!/bin/bash
#PBS -q jubio
#PBS -l nodes=10:ppn=24
#PBS -N dxtopo_3

cd /local/jubio/luisa/topologytracing/
#./assign_wo35 -pred_dist $OUTDIR/dist-pred-${BASENAME}.dat -beads_dist $OUTDIR/dist-calc-${BASENAME}.dat -startweights ${OUTDIR}/counts.txt -outdir ${OUTDIR}/assignmenttimetest > logassign
./assign_clean -pred_dist $OUTDIR/dist-pred-${BASENAME}.dat -beads_dist $OUTDIR/dist-calc-${BASENAME}.dat -startweights ${OUTDIR}/counts.txt -outdir ${OUTDIR}/assignmenttimetestw35tothe3 > logassign
EOF

chmod u+x jubio-assign-${BASENAME}.sh
qsub ./jubio-assign-${BASENAME}.sh

