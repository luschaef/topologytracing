#!/bin/bash

PDBFILE=$1
DATFILE=$2
#grep "^ATOM" ${PDBFILE} | awk '{ print $6, $7, $8 }' > ${DATFILE}
grep "^ATOM" ${PDBFILE} | awk '{ print $7, $8, $9 }' > ${DATFILE}
