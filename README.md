# topologytracing

Program to determine the topology of a protein given a cryo-EM density map and predicted inter-residue distances. Can also be used to assess topology of a given trace with regard to its agreement with predicted distances. 

Required input:  

                -density map in mrc format   
                
                -predicted inter-residue distances from TrRosetta[1] : <outdir>/seq.npz  (https://yanglab.nankai.edu.cn/trRosetta/)
                
                
Required software:  
                   
                   -Python 3.6 or newer  
                   
                   -LKH TSP solver[2] (http://webhotel4.ruc.dk/~keld/research/LKH/)  
                   
                   -DireX[3]  (http://schroderlab.org/software/direx/)  
                   
                   
                   DireX and LKH are needed only for building a trace, not for assessing. 
                   
                   C-scripts can be compiled with:
                   gcc assign.c -o assign -lm -openmp
                   gcc coords2pdb_n.c -o coords2pdb_n -lm 
                   
Usage:  
Building a trace:  

                        dxtopology.sh <MAP> <NRES> <BASENAME> <OUTDIR>  
                        
                        - MAP: mrc file of the density map  
                        
                        - NRES: Number of residues in the sequence of the protein of interest  
                        
                        - BASENAME: string to include in the name of all outputfiles  
                        
                        - OUTDIR: output directory, the seq.npz should be in this directory.   
                        
                       
Assessing a trace: 
       
                      dxtopologyV.sh <TRACE> <NRES> <BASENAME> <OUTDIR>
                      python assess_trace.py <BASENAME> <OUTDIR>
                      chimera <OUTDIR>/val.bild <OUTDIR>/val_alternatives.bild
                      
                      <TRACE> should be given as C_alpha-trace in PDB-format. In the .bild files C_alpha-atoms which might be placed unfavourable are 
                      depicted as gray spheres. Connections between beads which might be wrong are coloured red, proposed alternatives are depicted 
                      as green sticks.Expected positions of N- and C-terminus are blue and bright red.
                     
References:

                      [1] Yang, Jianyi, et al. "Improved protein structure prediction using predicted interresidue orientations." Proceedings of the National Academy of Sciences 117.3 (2020): 1496-1503.  
                      
                      [2] Helsgaun, Keld. "An effective implementation of the Lin–Kernighan traveling salesman heuristic." European journal of operational research 126.1 (2000): 106-130.   
                      
                      [3] G.F.Schröder, A.T.Brunger, and M.Levitt 'Combining Efficient Conformational Sampling with a Deformable Elastic Network Model Facilitates Structure Refinement at Low Resolution' Structure, Vol 15, 1630-1641 (2007)
