#!/bin/bash

OUTDIR="ribo-I-4A"
BASENAME="ribo-I-4A"
MAP="ribo-I-4A.mrc"
NRES="125"

#mkdir $OUTDIR
cp dxtraces-3.sh ${OUTDIR}/
chmod u+x ${OUTDIR}/dxtraces-3.sh
#initialize trace
cat << EOF > jubiojob-${BASENAME}.sh
#!/bin/bash
#PBS -q jubio
#PBS -l nodes=10:ppn=24
#PBS -N dxtopo_1

cd /local/jubio/luisa/topologytracing/$OUTDIR
./dxtraces-3.sh $MAP $NRES $OUTDIR > log
mv traces/refined-1.pdb ${BASENAME}-init-trace.pdb
EOF

chmod u+x jubiojob-${BASENAME}.sh
qsub ./jubiojob-${BASENAME}.sh

