#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 14:42:00 2022

@author: luisa
"""

import numpy as np
import sys 
import math
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

def calc_tracescores(tracelength,n, calcmat, predmat):
    print(tracelength)
    counts = np.zeros((n,n))
    partial_traces = [] 
    
    
    for startpos in range(n-tracelength+1):
        trace = [pos for pos in range(startpos, startpos+tracelength)]
        partial_traces.append(trace) 
        
    ########################################################### 
     ### Split predicted distance map in fitting submatrices ###
     ###########################################################
    sub_tr_length = tracelength
    inner_mats_pred = []
    outer_mats_pred =[]
    i=0
        
    while(i<=n-sub_tr_length):
         j = i+sub_tr_length      
         submat = pred_mat[i:j,:]
 
         inner_dists = submat[:,i:j]
         inner_dists = inner_dists[inner_dists >0]
     
         inner_mats_pred.append(inner_dists) 
         outer_dists= np.concatenate((submat[:,:i], submat[:,j:]),axis = 1)
  
         
         for k,row in enumerate(outer_dists):
             row = np.sort(row)
             outer_dists[k]= row
  
         outer_mats_pred.append(outer_dists) 
         i= i+1
         
         
    all_scores = []
    for t in partial_traces:
    
       
   
        
        #################################################
        ### split bead distance matrix in submatrices ###
        #################################################
        

        j = t[-1]+1
        i = t[0]
        submat = calc_mat[i:j,:]
 
        dists_inner = submat[:,i:j]
        dists_inner = dists_inner[dists_inner >0]
   
       
        dists_outer= np.concatenate((submat[:,:i], submat[:,j:]),axis = 1)

       
        for k,row in enumerate(dists_outer):
            row = np.sort(row)
            dists_outer[k]= row

       
         
      
          
      ################################################################################################       
      ### CaLCULATE SCORES AS CORRELATION BETWEEN calculated submatrices and predicted submatrices ###
      ################################################################################################
     
        scores = np.zeros(len(inner_mats_pred))
      
        
        for j in range((len(inner_mats_pred))):
        
          score_inner = np.sqrt(np.sum((dists_inner.ravel()-inner_mats_pred[j].ravel())**2))        
          score_inner_flip = np.sqrt(np.sum((np.flip(dists_inner).ravel()-inner_mats_pred[j].ravel())**2))
                 
          score_outer = np.sqrt(np.sum((dists_outer.ravel()-outer_mats_pred[j].ravel())**2))
          score_outer_flip = np.sqrt(np.sum((np.flipud(dists_outer).ravel()-outer_mats_pred[j].ravel())**2))
          if math.isnan(score_inner):
              scores[j] =  score_outer
          else:
              scores[j] =  min(((1/2*score_outer + 1/2*score_inner)),((1/2*score_outer_flip + 1/2*score_inner_flip)))
            
      
        all_scores.append(scores)
        
        
    all_scores = np.vstack(all_scores)
    all_scores = -1 * all_scores
   
 
    
    c_ind = int((tracelength-1)/2)
    c_beads = [trace[c_ind] for trace in partial_traces]
  
    maxs = []
    # keep best position for each bead
    maxscores = np.zeros_like(all_scores)
    for i,scores in enumerate(all_scores):
        l = list(scores)
        maxpos = l.index(max(scores)) + c_ind
        maxscores[i,maxpos-c_ind]=max(scores)
        maxs.append(maxpos)
    
    for i in range(len(maxs)):
        bead = c_beads[i]
        pos = maxs[i]
        counts[bead, pos] += 1

    return counts




basename = sys.argv[1]
outdir = sys.argv[2]

beads = np.loadtxt(f"{outdir}/{basename}-init-trace.dat")
#beads = np.loadtxt("/home/luisa/mnt/jubio3/distmaps/assignment_code/jan22_card/refined-1.dat")
n = len(beads)

pred_mat = np.loadtxt(f"{outdir}/dist-pred-{basename}.dat")
#pred_mat = np.loadtxt("/home/luisa/mnt/jubio3/distmaps/assignment_code/dist-pred-card.dat")
pred_mat = pred_mat.reshape(n,n)
calc_mat = np.loadtxt(f"{outdir}/dist-calc-{basename}.dat")
#calc_mat = np.loadtxt("/home/luisa/mnt/jubio3/distmaps/assignment_code/jan22_card/card-dist-calc.dat")
calc_mat = calc_mat.reshape((n,n))     



partial_traces = [] 

sum_counts = np.zeros((n,n))
maxlength = int(n-0.1*n)
if (maxlength % 2 == 0):
    maxlength = maxlength +1


 
results = Parallel(n_jobs=4)(delayed(calc_tracescores)(tracelength,n,calc_mat,pred_mat) for tracelength in range(maxlength,3,-2))    
for result in results:
    sum_counts = sum_counts+result    
 
np.savetxt(f"{outdir}/counts.txt", sum_counts.ravel(), fmt="%d")

plt.figure()
plt.imshow(sum_counts)
    
    
    
 
