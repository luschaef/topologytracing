#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 09:20:41 2022

@author: luisa
"""

import numpy as np
import sys

outdir = sys.argv[2]
basename = sys.argv[1]
indmax = np.loadtxt(f"{outdir}/assignment/ind_max_999.dat",dtype=int)
beads = np.loadtxt(f"{outdir}/{basename}-init-trace.dat")
outbeads = np.array([beads[i] for i in indmax])
dists = [np.sqrt(np.sum((outbeads[i] -outbeads[i-1])**2)) for i in range(1,len(outbeads))]

outbeads = np.array([beads[i]for i in indmax])
np.savetxt(f"{outdir}/out-raw-{basename}.dat", outbeads, fmt="%.3f")
n_raw = len(outbeads)
f = open(f"{outdir}/n_raw.txt", "w")
f.write(str(n_raw))
f.close()


#remove outliers
clean_inds = [indmax[i] for i in range(0,len(dists)-1) if (dists[i] < 10 and dists[i+1] < 10) ]
clean_outbeads = np.array([beads[i]for i in clean_inds])
clean_dists = [np.sqrt(np.sum((clean_outbeads[i] - clean_outbeads[i-1])**2)) for i in range(1,len(clean_outbeads))]
np.savetxt(f"{outdir}/out-cleaned-{basename}.dat", clean_outbeads, fmt="%.3f")
n_clean = len(clean_outbeads)
f = open(f"{outdir}/n_clean.txt", "w")
f.write(str(n_clean))
f.close()


#"smoothing step 1: extract diagonals, partial traces" 
idists = np.abs(np.diff(np.array(clean_inds)))
ithr = 10 #je größer desto mehr schwankung ist innerhalb einer diagonalen erlaubt 
pt_1 = []
trace = []
for j in range(len(idists)):
    if idists[j] <= ithr:
        trace.append(clean_inds[j])
       
    else:
        trace.append(clean_inds[j])
        pt_1.append(trace)
        trace = []
       
if np.abs(clean_inds[-1] - clean_inds[-2]) <=ithr:
    trace.append(clean_inds[-1])        
pt_1.append(trace)

## dont regard single beads
pt_2 = []        
for trace in pt_1:
    if len(trace)>1:
        pt_2.append(trace)

fulltrace1 = []
for trace in pt_2:
    fulltrace1 = fulltrace1 + list(trace)

smooth_outbeads = [beads[i] for i in fulltrace1]
np.savetxt(f"{outdir}/out-smoothed-{basename}.dat", smooth_outbeads, fmt="%.3f")

n_smooth = len(smooth_outbeads)
f = open(f"{outdir}/n_smooth.txt", "w")
f.write(str(n_smooth))
f.close()


#sort traces
pt_3 = []
for trace in pt_2:    
    if np.mean(np.diff(trace)) > 0:
        sorted_trace = np.sort(trace)
    else:
        sorted_trace = np.flip(np.sort(trace))
    pt_3.append(sorted_trace)

fulltrace2 = []
for trace in pt_3:
    fulltrace2 = fulltrace2 + list(trace)    
sorted_outbeads = [beads[i] for i in fulltrace2]
np.savetxt(f"{outdir}/out-sorted-{basename}.dat", sorted_outbeads, fmt="%.3f")

n_sorted = len(sorted_outbeads)
f = open(f"{outdir}/n_sorted.txt", "w")
f.write(str(n_sorted))
f.close()

#fill traces
notassigned = [ind for ind in indmax if ind not in fulltrace2]
pt_4 = pt_3
for ind in notassigned:
    for i,trace in enumerate(pt_3):
        if (min(trace) < ind < max(trace)):
            #print(ind, trace)
            trace = np.append(trace,ind)
            notassigned.remove(ind)
            #print(trace)
            if np.mean(np.diff(trace)) > 0:
                sorted_trace = np.sort(trace)
            else:
                sorted_trace = np.flip(np.sort(trace))
            pt_4[i] = sorted_trace  
            break
        
fulltrace3 = []
for trace in pt_4:
    fulltrace3 = fulltrace3 + list(trace)    
filled_outbeads = [beads[i] for i in fulltrace3]
np.savetxt(f"{outdir}/out-filled-{basename}.dat", filled_outbeads, fmt="%.3f")

n_filled = len(filled_outbeads)
f = open(f"{outdir}/n_filled.txt", "w")
f.write(str(n_filled))
f.close()
