#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import sys

basename=sys.argv[1]
outdir=sys.argv[2]



f = open(f"{outdir}/val.bild", "w")
g = open(f"{outdir}/val_alternatives.bild", "w")
origbeads = np.loadtxt(f"{outdir}/{basename}-init-trace.dat")
filledbeads = np.loadtxt(f"{outdir}/out-filled-{basename}.dat")
n=len(origbeads)
js = []
for i, bead in enumerate(filledbeads):
    j = np.nonzero((origbeads == bead))[0]
    #print(j)
    if len(j > 3):
        for k in j:
            c = (j == k).sum()
            if c == 3:
                l = k
    else:
        l= j[0]  
    js.append(l)
    

missing = []
wrongs = []
wrongs_m = []
inds = np.arange(n)
for ind in inds:
    if ind > n-1:
        continue;
    next_ind = ind+1
    j = np.nonzero(js == ind)[0]
    if len(j)>0:
        j = j[0]
        if j > len(js)-2:
            continue;
        if js[j+1] != next_ind and js[j-1] != next_ind and len(np.nonzero(js == next_ind)[0]) > 0 :
             wrongs.append((ind, next_ind))
    else:
        missing.append(ind)
        if ind < n-2:
            wrongs_m.append((ind, ind+1))
        if ind > 0:
            wrongs_m.append((ind-1,ind))

        
wrongs = set(wrongs)
wrongs_m = set(wrongs_m)

###find alternative:
alts= []
for wrong in wrongs:
    ind = wrong[0]
    ind_spot = np.nonzero(js == ind)[0][0]
    if ind_spot < (len(js)-1):
        if abs(js[ind_spot+1] -ind) > 1 :
            alts.append((ind, js[ind_spot+1]))
    elif ind_spot > 0:
        alts.append((js[ind_spot-1],ind))

    ind = wrong[1]
    ind_spot = np.nonzero(js == ind)[0][0]
    if ind_spot < (len(js)-1):
        if abs(js[ind_spot+1] -ind) > 1 :
            alts.append((ind, js[ind_spot+1]))
    elif ind_spot > 0:
        alts.append((js[ind_spot-1],ind))

alts = set(alts)
   
    
alts_w = []
for w in wrongs_m:
   if w[0] in missing and w[1] in missing:
       continue
   for ind in w:
       if ind not in missing:
           nm = ind
           print(nm,w)
   ind_nm = np.nonzero(js == nm)[0][0]
   if ind_nm < (len(js)-1):
       if abs(js[ind_nm+1] -nm) > 1 :
           alts_w.append((nm, js[ind_nm+1]))
   elif ind_nm > 0:
       alts_w.append((js[ind_nm-1], nm))
       
alts_w = set(alts_w)   
        
   
#write bild file 


f.write(".color gray \n")

for m in missing:
    m = m
    x = origbeads[m,0]
    y = origbeads[m,1]
    z = origbeads[m,2]
    f.write(".sphere "+str(x)+" " +str(y)+" " +str(z)+" 0.8 \n") 

f.write(".color 0.67 0.29 0.29 \n") 
for w in wrongs:
    p1 = w[0] 
    p2 = w[1] 
    x1 = origbeads[p1,0]
    y1 = origbeads[p1,1]
    z1 = origbeads[p1,2]
    x2 = origbeads[p2,0]
    y2 = origbeads[p2,1]
    z2 = origbeads[p2,2]
    
    f.write(".cylinder "+str(x1)+" " +str(y1)+" " +str(z1)+ " " + str(x2)+" " +str(y2)+" " +str(z2)+ " 0.2 \n")
    f.write(".sphere "+str(x1)+" " +str(y1)+" " +str(z1)+" 0.5 \n") 
    f.write(".sphere "+str(x2)+" " +str(y2)+" " +str(z2)+" 0.5 \n") 

f.write(".color 0.19 0.34 0.35 \n")
for w in wrongs_m:
    p1 = w[0] 
    p2 = w[1] 
    x1 = origbeads[p1,0]
    y1 = origbeads[p1,1]
    z1 = origbeads[p1,2]
    x2 = origbeads[p2,0]
    y2 = origbeads[p2,1]
    z2 = origbeads[p2,2]
    
    f.write(".cylinder "+str(x1)+" " +str(y1)+" " +str(z1)+ " " + str(x2)+" " +str(y2)+" " +str(z2)+ " 0.2 \n")
    f.write(".sphere "+str(x1)+" " +str(y1)+" " +str(z1)+" 0.5 \n") 
    f.write(".sphere "+str(x2)+" " +str(y2)+" " +str(z2)+" 0.5 \n")    

    
n_term = js[0] 
x = origbeads[n_term,0]
y = origbeads[n_term,1]
z = origbeads[n_term,2]
f.write(".color medium blue \n")
f.write(".sphere "+str(x)+" " +str(y)+" " +str(z)+" 0.8 \n") 
  

c_term = js[-1] 
x = origbeads[c_term,0]
y = origbeads[c_term,1]
z = origbeads[c_term,2]
f.write(".color red \n")
f.write(".sphere "+str(x)+" " +str(y)+" " +str(z)+" 0.8 \n")   
f.close()


g.write(".color 0.34 0.65 0.45 \n") 
for a in alts:
    p1 = a[0] 
    p2 = a[1] 
    x1 = origbeads[p1,0]
    y1 = origbeads[p1,1]
    z1 = origbeads[p1,2]
    x2 = origbeads[p2,0]
    y2 = origbeads[p2,1]
    z2 = origbeads[p2,2]
    
    g.write(".cylinder "+str(x1)+" " +str(y1)+" " +str(z1)+ " " + str(x2)+" " +str(y2)+" " +str(z2)+ " 0.2 \n")
    g.write(".sphere "+str(x1)+" " +str(y1)+" " +str(z1)+" 0.5 \n") 
    g.write(".sphere "+str(x2)+" " +str(y2)+" " +str(z2)+" 0.5 \n") 
    
g.write(".color 0.34 0.65 0.45 \n")    
for a in alts_w:
    p1 = a[0] 
    p2 = a[1] 
    x1 = origbeads[p1,0]
    y1 = origbeads[p1,1]
    z1 = origbeads[p1,2]
    x2 = origbeads[p2,0]
    y2 = origbeads[p2,1]
    z2 = origbeads[p2,2]
    
    g.write(".cylinder "+str(x1)+" " +str(y1)+" " +str(z1)+ " " + str(x2)+" " +str(y2)+" " +str(z2)+ " 0.2 \n")
    g.write(".sphere "+str(x1)+" " +str(y1)+" " +str(z1)+" 0.5 \n") 
    g.write(".sphere "+str(x2)+" " +str(y2)+" " +str(z2)+" 0.5 \n") 
    
g.close()
