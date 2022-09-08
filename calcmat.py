#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

outdir = sys.argv[2]
basename = sys.argv[1]
beads = np.loadtxt(f"{outdir}/{basename}-init-trace.dat")
n = len(beads)
pred_mat = np.loadtxt(f"{outdir}/dist-pred-{basename}.dat")

pred_mat = pred_mat.reshape(n,n)
calc_mat = np.zeros_like(pred_mat)
for k,bead in enumerate(beads):
    for l,other_bead in enumerate(beads):
        calc_mat[k,l] = np.sqrt(np.sum((bead-other_bead)**2))
       
sort_inds = np.argsort(pred_mat.ravel())
sort_inds2 = np.argsort(calc_mat.ravel())

calc_mat.ravel()[sort_inds2] = pred_mat.ravel()[sort_inds] 

np.savetxt(f"{outdir}/dist-calc-{basename}.dat", calc_mat.ravel(), fmt="%d") 
