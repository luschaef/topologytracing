#!/bin/bash


## PROGRAM TO GENERATE A CA-TRACE FROM A
## DENSITY MAP
## 2016-08-14

## BY GUNNAR F SCHROEDER
##
## needs:
## EMAN2
## DireX
## Chimera

outdir=traces
inputmap=$1
# NUMBER OF RESIDUES
NRES=$2

## threshold of locally normalized map "in.mrc"ls
## used for e2pathwalker.py
MAP_THRESHOLD=1.0

# NUMBER OF TRACES USED TO BUILD CONNECTION HISTOGRAMS
# THE MORE RESIDUES THE LARGER SHOULD THIS VALUE BE.
NAVG=10
# NUMBER OF TOTAL TRACES TO GENERATE
NTRACES=1

# DEFINE RESOLUTION RANGES FOR DIREX REFINEMENT
CUR_MAP_DMIN=5
MAP_CV_DMIN=4





#+++++++++++++++++++++++++++++++++++++++++++++++++
# THINGS BELOW THIS LINE DO NOT NEED TO BE CHANGED
#+++++++++++++++++++++++++++++++++++++++++++++++++

### 1=YES, 0=NO
PREPARE_MAP=1
RUN_PATHWALKER=1
REFINE_TRACE=1

# We typically build twice the number of beads
# as there are residues.
NRES_FAC=1
NRES=$(( NRES * NRES_FAC ))

mkdir ${outdir}
rm -rf  ./.tmp-trace-*.pdb*
##################
### PREPARE MAP###
##################

if (( ${PREPARE_MAP} ))
then
echo "Prepare Map..."

#shift map to origin 
cat << EOF > move-map.com
open ${inputmap}
volume #0 originIndex 0,0,0
volume #0 save in0.mrc
EOF
chimera --nogui --silent move-map.com
rm move-map.com
e2proc3d.py in0.mrc in1.mrc --process=normalize

e2proc3d.py in1.mrc in.mrc --process=normalize.local:radius=10:threshold=3.0 --process=filter.lowpass.tophat:cutoff_freq=0.2
e2proc3d.py in.mrc proc.mrc --process=threshold.belowtozero 
#rm in1.mrc

cat << EOF > prepare-map.com
open proc.mrc
vop gaussian #0 sDev 1.0
vop ridges #1 level 0.5
vop gaussian #2 sDev 1.0
volume #3 save out.mrc
EOF
#cat << EOF > prepare-map.com
#open proc.mrc
#vop ridges #1 level 0.5
#volume #3 save out.mrc
#EOF
chimera --nogui --silent prepare-map.com
rm prepare-map.com
#rm out.mrc
e2proc3d.py out.mrc out1.mrc --process=threshold.belowtozero:minval=1.0
fi


######################
### RUN PATHWALKER ###
######################

if (( ${RUN_PATHWALKER} ))
then
for n in `seq 1 ${NTRACES}`
do
echo "Working on Trace $n"
### MAKE BEAD MODEL
if (( 1 ))
then
echo "Make Bead Model..."
dxbeadgen out1.mrc beadcoordinate_val.pdb ${NRES} 3 0 
#rm out1.mrc

 
## REFINE BEAD MODEL
cat << EOF > dx.par
nsteps     = 100
annealing     = 100
concoord_damp = 1.000000
nviol      = 100000
pert_fac = 0.10000
use_den = no
map_strength = 0.0010000
map_damp_restraints = 100
cur_map_dmin = ${CUR_MAP_DMIN}
map_cv_dmin = ${MAP_CV_DMIN}
min_cycles = 100
no_bonds = yes
repel_damp = 100
#repel_all_atoms = 0.000
repel_all_atoms = 0.005
vdw_lb_dif = 1.0
vdw_ub_dif = 1.0
use_planar_constraints = no
use_improper_constraints = no
scale_difmap = 1.0
bead_mode = yes
EOF
direx -pdb beadcoordinate_val.pdb -f dx.par -o c.pdb -map proc.mrc

fi
#rm dx.par


for i in `seq 1 ${NAVG}`
do
### GENERATE TRACE
echo "running pathwalker"
#/home/luisa/e2pathwalker_gunnar.py --output=.tmp-trace-${i}.pdb --solver=lkh  --mapfile=in1.mrc --mapthresh=${MAP_THRESHOLD} --mapweight=0.6  --average=0 --noise=0.1 --subunit=1 --overwrite --dmin=0.0 --dmax=8.0  c.pdb
e2pathwalker_customized.py --output=.tmp-trace-${i}.pdb --solver=lkh  --mapfile=${inputmap} --mapthresh=${MAP_THRESHOLD} --mapweight=0.6  --average=0 --noise=0.1 --subunit=1 --overwrite --dmin=0.0 --dmax=8.0  c.pdb
#mv .tmp-trace-${i}.pdb trace-${i}.pdb 
done


#rm c.pdb

#######
#./mkhist.sh .tmp-trace-*.pdb
rm -rf .tmp.hist

N=`ls -1 .tmp-trace-*.pdb | wc | awk '{print $1}'`

if (($N<1))
then
  echo "no input files found for mkhist!"
  exit
fi

COUNT=0
for i in .tmp-trace-*.pdb
do
awk ' BEGIN {k=0} /^ATOM/ {if(k>0){ if(old<$6){print old+1,$6+1;} else {print $6+1,old+1}}; old=$6; k++} ' ${i} >> .tmp.hist
NATOMS=`grep "^ATOM" ${i} | wc | awk '{print $1}'`
COUNT=$((COUNT+1))
done

echo "COUNT= $COUNT"
echo "NATOMS= $NATOMS"

# WRITE ALL PAIRS
for ((i=1; i<$NATOMS+1; i++))
do
  for ((j=i+1; j<=$NATOMS+1; j++))
  do
     echo $i $j >> .tmp.hist
  done
done

# WRITE HEADER
cat <<EOF > problem.tsp
NAME: this
TYPE: TSP
COMMENT: this
DIMENSION: $((NATOMS+1))
EDGE_WEIGHT_TYPE: EXPLICIT
EDGE_WEIGHT_FORMAT: UPPER_ROW
FIXED_EDGES_SECTION
-1
EDGE_WEIGHT_SECTION
EOF

sort -n -k1,1 -k2,2 .tmp.hist  | uniq -c | awk '{print $2,$3,$1}' | awk '{if (($1<='${NATOMS}') && ($2<='${NATOMS}')) {print int('$COUNT'*100./$3)} else {print "0"} }'  >> problem.tsp
rm .tmp.hist

cat <<EOF > params.lkh
PROBLEM_FILE = problem.tsp
OUTPUT_TOUR_FILE = lkh.out
RUNS = 10
TRACE_LEVEL = 5
EOF

### RUN LKH
LKH params.lkh

##############

rm -rf ${outdir}/trace-${n}.pdb

NATOMS=`grep " CA " .tmp-trace-1.pdb | wc | awk '{print $1}'`

awk 'NR>6 && NR<'${NATOMS}'+1+6+1 {print}' lkh.out > .tmp.list
awk 'NR>6 && NR<'${NATOMS}'+1+6+1 {print}' lkh.out >> .tmp.list

awk ' BEGIN {s=0} { if (s==1) {print}; if ($1=="'$((NATOMS+1))'") {s=s+1} } ' .tmp.list > .tmp.list2

for i in `cat .tmp.list2`
do
   if (( $i < $NATOMS + 1))
   then
     i=$((i-1))
     awk '/^ATOM/ && $6=="'${i}'"' .tmp-trace-1.pdb >> ${outdir}/trace-${n}.pdb
   fi
done

echo "END" >> ${outdir}/trace-${n}.pdb

#rm .tmp.list
#rm .tmp.list2


awk '/^ATOM/ { printf "ATOM  %5i%-11s%4i%s\n",NR,substr($0,12,11),NR,substr($0,27,90);} ' ${outdir}/trace-${n}.pdb > .tmp.pdb
echo "END" >> .tmp.pdb
mv .tmp.pdb ${outdir}/trace-${n}.pdb


done

fi  #zu Pathwalker


####################
### REFINE TRACE ###
####################


if (( ${REFINE_TRACE} ))
then
echo "Refine Traces..."
for n in `seq 1 ${NTRACES}`
do

echo "Refining Trace ${n}"
### GENERATE RESTRAINTS
awk 'BEGIN {n=0} /^ATOM/ {if (n % '${NRES_FAC}' ==0) {print}; n=n+1 }'  ${outdir}/trace-${n}.pdb > tmp.pdb
awk '/^ATOM/ { printf "ATOM  %5i%-11s%4i%s\n",NR,substr($0,12,11),NR,substr($0,27,90);} ' tmp.pdb > tmp2.pdb
echo "END" >> tmp2.pdb

AVG_BOND_LENGTH=`awk 'BEGIN {print 3.8/'${NRES_FAC}' }'`
echo "AVG_BOND_LENGTH= ${AVG_BOND_LENGTH}"
## BONDS
awk ' BEGIN {a=0;n=0;} /^ATOM/ {if (n>0) { print $2, a," 3.8  1.00  1.00" }; a=$2;n=n+1 }' tmp2.pdb > .tmp.restraints
## ANGLES  ( 5.2 - 7.2 A )
awk ' BEGIN {a1=0;a2=0;n=0;} /^ATOM/ {if (n>1) { print $2, a2," 6.0  0.10  1.00" }; a2=a1;a1=$2;n=n+1 }' tmp2.pdb >> .tmp.restraints
wc .tmp.restraints | awk '{print $1}' > ${outdir}/restraints-${n}.dat
cat .tmp.restraints >> ${outdir}/restraints-${n}.dat
rm .tmp.restraints

### REFINE TRACE
cat << EOF > dx-refine-trace.par
nsteps     = 100
annealing     = 100
concoord_damp = 1.000000
nviol      = 0
pert_fac = 0.10000
use_den = no
map_strength = 0.0010000
map_damp_restraints = 100
cur_map_dmin = ${CUR_MAP_DMIN}
map_cv_dmin = ${MAP_CV_DMIN}
min_cycles = 100
no_bonds = yes
repel_damp = 100
repel_all_atoms = 0.000
#repel_all_atoms = 0.0005
vdw_lb_dif = 1.0
vdw_ub_dif = 1.0
use_planar_constraints = no
use_improper_constraints = no
scale_difmap = 1.0
bead_mode = yes
EOF

/home/gschroed/bin/direx -pdb tmp2.pdb -f dx-refine-trace.par -o ${outdir}/refined-${n}.pdb -map in.mrc -expd ${outdir}/restraints-${n}.dat 

### generate inverse trace
#grep ^ATOM ${outdir}/refined-${n}.pdb | tac > ${outdir}/refined-${n}-reverse.pdb
#echo "END" >>  ${outdir}/refined-${n}-reverse.pdb

done

fi

rm restraints*
rm ./tmp.pdb ./tmp2.pdb
rm ./.tmp-trace-*.pdb*
rm ./problem.tsp
rm ./params.lkh
rm ./lkh.out
#rm out1.mrc 
rm dx-refine-trace.par
