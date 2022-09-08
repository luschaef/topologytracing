import numpy as np
import matplotlib.pyplot as plt
import sys

filename = sys.argv[1]
basename= sys.argv[2] 
outdir = sys.argv[3]
n = int(sys.argv[4])

data = np.load(filename, allow_pickle=True)
lst = data.files


dist = data['dist']
a=np.zeros((n,n))

for i in range(n):
    for j in range(n):
        if (i==j):
            a[i,j]= 0
        else:
            a[i,j]= np.argmax(dist[i,j,1:])



np.savetxt(f'{outdir}/dist-pred-{basename}.dat', a.ravel(), fmt='%d') 


