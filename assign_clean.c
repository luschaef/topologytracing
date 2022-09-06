#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fftw3.h>
#include <time.h>

#define MAX_STRING_SIZE 9999
#define DIST_CUTOFF 15
#define MAX_NEIGHBORS 20
#define N_SMOOTH_ITERATIONS 0


double kappa_dist_weight = 1.00;
double kappa_dist_scale = 1.;
char fname[MAX_STRING_SIZE];
char beadsfname[MAX_STRING_SIZE];
char preddistfname[MAX_STRING_SIZE];
char outdir[MAX_STRING_SIZE];
char outfname[MAX_STRING_SIZE];
char distfname[MAX_STRING_SIZE];
char startweightsfname[MAX_STRING_SIZE];
int NSTEPS=1000;
int arg_count;
int n, nseq;
int i,j,k,l,h,m,g;
int kk;
int count;
double nsqr;


double grad_factor = 1000.0;
double max_gradient = 1.0;

double score;
double *s;
double f1,f2,f3;
double max_value, min_value;
FILE *fp;

double *dist_matrix_pred;
double *dist_matrix;
double *tmp_dist_matrix;
double *w, *w_copy, *grad_w, *fac_grad;
double *p;



int *ind_max;
int ind_max_n, ind_max_nseq;



static inline double sqr( double x )
{
  return (x * x);
}

// ************************************ MAIN *******************************************************************************************

int main (int argc, char *argv[])
{
clock_t begin = clock();
//=============================
// PARSE COMMAND LINE ARGUMENTS
// ============================

  if ( argc < 4 ) printf("Too few arguments");
  arg_count = 1;
  while ( arg_count < argc ){
  //printf("Argci: %d\n", argc);
  //printf("arg_count: %d\n", arg_count);

   
     if ( ! strcmp(argv[arg_count],"-beads_dist") )
     {
         strcpy(distfname,argv[arg_count+1]);
         printf("Bead Distances .............................. : %s\n",distfname); 
         arg_count+=2;

     } else if ( ! strcmp(argv[arg_count],"-pred_dist") )
     {
         strcpy(preddistfname,argv[arg_count+1]);
         printf("Predicted distances .......................... : %s\n",preddistfname); 
         arg_count+=2;
	
     } else if ( ! strcmp(argv[arg_count],"-startweights") )
     {
         strcpy(startweightsfname,argv[arg_count+1]);
         printf("Initial weights .......................... : %s\n",startweightsfname); 
         arg_count+=2;

     } else if ( ! strcmp(argv[arg_count], "-outdir") )
     {
         strcpy(outdir,argv[arg_count+1]);
         printf("Output directory ...........................: %s\n", outdir);
         mkdir(outdir, S_IRWXU | S_IRWXG );
         arg_count+=2;}
     
      	
}

// ===================================
// READ DISTANCE MATRIX FROM trRosetta
// ===================================

printf("Reading predicted Distances\n");

     if ( (dist_matrix_pred = (double *) calloc(1000000,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
     }

     strcpy(fname,preddistfname);
     if ( (fp = fopen(fname,"r")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }

     i=0;
     while ( ! feof(fp) ){
	fscanf(fp, "%lf", &dist_matrix_pred[i]);
	i++;
}

     fclose(fp);
     nsqr = (double) i ;

     nseq =  (int) 	sqrt(nsqr);
     printf("Done,nseq=%d\n", nseq);

// ===================================
// READ DISTANCE MATRIX FROM BEADS
// ===================================

printf("Reading calculated Distances\n");

     if ( (dist_matrix = (double *) calloc(1000000,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
     }

     strcpy(fname,distfname);
     if ( (fp = fopen(fname,"r")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }

     i=0;
     while ( ! feof(fp) ){
	fscanf(fp, "%lf", &dist_matrix[i]);
	i++;
}

     fclose(fp);
     n = (double) i ;

     n =  (int) 	sqrt(nsqr);
     printf("Done,n=%d\n", nseq);


//==============================
// ASSIGNMENT WEIGHTS
// ============================
fprintf(stderr,"Initialize weight matrix\n");
    if ( (w= (double *) calloc(nseq*n,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
    }
    if ( (w_copy= (double *) calloc(nseq*n,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
    }

     strcpy(fname,startweightsfname);
     if ( (fp = fopen(fname,"r")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }

     i=0;
     while ( ! feof(fp) ){
	fscanf(fp, "%lf", &w[i]);
	i++;
}
    //for (i=0;i<nseq*n;i++) w[i]*= 10;
    for (i=0;i<nseq*n;i++) w[i]+= 1.0/n;


    if ( (grad_w= (double *) calloc(nseq*n,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
    }
// NORMALIZE WEIGHTS
// w is now squared
  for (i=0;i<nseq;i++) {
    // CALC SUM
    f1= 0.0;
    for (j=0;j<n;j++)  {
       f1+= sqr(w[j*nseq+i]);
    }
    // NORMALIZE
    if (f1 > 0) f1=1.0/f1;
    else f1=0;
    //f1=1.0/f1;
    for (j=0;j<n;j++)  {
       w[j*nseq+i]*= sqrt(f1);
    }
  } 


if (nseq>n) { 
   printf("nseq>n, this should not happen. You need more beads.\n"); 
   exit(1);
}

// CREATE ind_max array
if ( (ind_max= (int *) calloc(nseq,sizeof(int)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
}


if ( (s = (double *) calloc(NSTEPS,sizeof(double)) ) ==  NULL ) {
        fprintf(stderr,"ERROR: could not allocate memory \n");
        exit(1);
}



// ===============
//    MAIN LOOP 
// ===============

for (kk=0;kk<NSTEPS;kk++) {

fprintf(stderr,"step %i\n",kk);

  // ===============
  // CALCULATE SCORE
  // ===============

  fprintf(stderr,"calc score...\n");
  // CALCULATE DISTANCE MATRIX RMS DIFFERENCE
  score= 0.0;

#pragma omp parallel private(i,j,k,l,f1,f2) reduction(+:score)
{

#pragma omp for
  for (i=0;i<nseq;i++) {
    for (j=0;j<nseq;j++) {
      for (k=0;k<n;k++) {
        for (l=0;l<n;l++) {
		
             	 f1 =0.;
	
                 f2 =  (kappa_dist_weight * sqr( dist_matrix[k*n+l] - dist_matrix_pred[i*nseq+j])) + 1.;
                 f1 = sqr(w[l*nseq+j]) * sqr(w[k*nseq+i]) *  kappa_dist_weight /  f2  ;
		
		 score -=f1;
        }
      }
    }
  }

} // zu openmp


  fprintf(stderr,"Score: %f\n",score/((double) n*nseq));
  s[kk]= score/((double) n*nseq);

  fprintf(stderr,"Calc gradient...\n");
  // ========================================
  // CALCULATE GRADIENT ON ASSIGNMENT WEIGHTS
  // ========================================
  for (h=0;h<nseq;h++) {
    for (m=0;m<n;m++) {
        grad_w[m*nseq+h]= 0.0;

        f1=0.0;
        #pragma omp parallel private(i,k,f2) reduction(+:f1)
        {
            #pragma omp for
            for (i=0;i<nseq;i++) {
              for (k=0;k<n;k++) {
                       f2= (kappa_dist_scale * sqr( dist_matrix[k*n+m] - dist_matrix_pred[i*nseq+h] )) + 1.;
                      f1 +=  sqr(w[k*nseq+i]) * kappa_dist_weight / f2 ;
              }
            }
        } // zu openmp
  
        //grad_w[m*nseq+h] = -2. * 2. * f1 * w[m*nseq+h] + 2. * w[m*nseq+h] ;
        grad_w[m*nseq+h] = -2. * 2. * f1 * w[m*nseq+h] + 2. * w[m*nseq+h]*w[m*nseq+h]*w[m*nseq+h] ;
    }
  }


  fprintf(stderr,"Scale gradient...\n");
  // ==============
  // SCALE GRADIENT
  // ==============
  
  // PRE-SCALE GRADIENT
  for (i=0;i<nseq;i++)  {
    for (k=0;k<n;k++)  {
       grad_w[k*nseq+i] /= (double) n * nseq;
    }
  }

  // ==============
  // UPDATE WEIGHTS
  // ==============
  for (i=0;i<nseq*n;i++) {
      w[i]-= grad_factor * grad_w[i];
	f2 = w[i];
  }
  





  fprintf(stderr,"Normalize weights...\n");
  // NORMALIZE WEIGHTS
  // in n-direction
  for (i=0;i<nseq;i++) {
    // CALC SUM
    f1= 0.0;
    for (j=0;j<n;j++)  {
       f1+= sqr(w[i*nseq+j]);
    }
    // NORMALIZE
    if (f1 > 0) f1=1.0/f1;                                                                                                                                                       
    else f1=0;  
//    f1=sqrt(f1);
//    f1=1.0/f1;
    for (j=0;j<n;j++)  {
       w[i*nseq+j]*= sqrt(f1);
    }
  } 
//
  // in nseq-direction
  for (i=0;i<nseq;i++) {
    // CALC SUM
    f1= 0.0;
    for (j=0;j<n;j++)  {
       f1+= sqr(w[j*nseq+i]);
    }
    // NORMALIZE
    if (f1 > 0) f1=1.0/f1;                                                                                                                                                       
    else f1=0;  
    //f1=sqrt(f1);
    //f1=1.0/f1;
    for (j=0;j<n;j++)  {
       w[j*nseq+i]*= sqrt(f1);
    }
  } 



  // just print
      printf("w  grad_w (0,0) =  %f   %f\n", sqr(w[0*nseq+0]), grad_w[0*nseq+0]);
      printf("w  grad_w (1,1) =  %f   %f\n", sqr(w[1*nseq+1]), grad_w[1*nseq+1]);
      printf("w  grad_w (0,1) =  %f   %f\n", sqr(w[0*nseq+1]), grad_w[0*nseq+1]);
      printf("w  grad_w (1,0) =  %f   %f\n", sqr(w[1*nseq+0]), grad_w[1*nseq+0]);

    f1= 0.0;
    for (j=0;j<n;j++)  f1+= sqr(w[j*nseq+0]);
    printf("sum(w(0))   %f\n",f1);

    f1= 0.0;
    for (j=0;j<n;j++)  f1+= sqr(w[j*nseq+1]);
    printf("sum(w(1))   %f\n",f1);

    f1= 0.0;
    for (j=0;j<n;j++)  f1+= sqr(w[j*nseq+20]);
    printf("sum(w(20))   %f\n",f1);



// ================================
// FIND ASSIGNMENT FROM MAX WEIGHTS
// ================================
// REMINDER: ind_max input is nseq and output is n

  fprintf(stderr,"Find Assignment..\n");
  // MAKE WEIGHT COPY
  for (i=0;i<n*nseq;i++) w_copy[i] = w[i];
  // or faster with memcpy
  // memcpy(w_copy, w, n*nseq*sizeof(double));
  
  for (i=0;i<nseq;i++) {
     max_value= -9999;
     // SEARCH MAX VALUE IN ENTIRE MATRIX
     for (k=0;k<n;k++) {
       for (j=0;j<nseq;j++) {
             if ( max_value < w_copy[k*nseq+j] ) {
                   max_value=  w_copy[k*nseq+j]; 
                   ind_max_n= k;
                   ind_max_nseq= j;
             }
       }
     }
     // SET ind_max
     ind_max[ind_max_nseq] = ind_max_n;

     // SET ENTIRE ROW AND COLUMN TO ZERO   
     // TO ENSURE UNIQUE ASSIGNMENT
     for (k=0;k<n;k++)  w_copy[k*nseq + ind_max_nseq] = -1.0;
     for (j=0;j<nseq;j++)  w_copy[ind_max_n*nseq + j] = -1.0;
  }   


// =======================
//    OUTPUT
// =======================

  fprintf(stderr,"Write output...\n");
if (0 == 0) { 
// write ind_max to file for checking
   sprintf(fname,"%s/ind_max_%d.dat", outdir,kk);
   //sprintf(fname,"ind_max_%05i.dat",kk);
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
   for (i=0;i<nseq;i++) {
      fprintf(fp,"%i\n",ind_max[i]);
   }
   fclose(fp);

// WRITE NEW DISTANCE MATRIX

   sprintf(fname,"%s/d.dat",outdir);
   //strcpy(fname,"d.dat");
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<nseq;i++) 
       for (j=0;j<nseq;j++) 
           fprintf(fp,"%f\n",dist_matrix[ind_max[i]*nseq + ind_max[j]]);
     fclose(fp);



// WRITE WEIGHTS
//   sprintf(fname,"%s/weights_%d.dat", outdir, kk);
   sprintf(fname,"%s/weights_%d.dat", outdir,kk);
   if ( (fp = fopen(fname,"w")) == NULL)
   {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
   }
// Achtung so ist w später transformiert zu counts und so! besser for j, for i
   for (i=0;i<nseq;i++) 
     for (j=0;j<n;j++) 
           fprintf(fp,"%f\n",w[j*nseq + i]);
   fclose(fp);

// WRITE GRADIENT
   //sprintf(fname,"%s/gradient_%d.dat", outdir,kk);
   sprintf(fname,"%s/gradient_%d.dat", outdir,kk);
   if ( (fp = fopen(fname,"w")) == NULL)
   {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
   }
   for (i=0;i<nseq;i++) 
     for (j=0;j<n;j++) 
           fprintf(fp,"%f\n",grad_factor *grad_w[j*nseq + i]);
   fclose(fp);

}




// ==================
//  END OF MAIN LOOP
// =================
}  // zu kk NSTEPS
	
// WRITE INPUT DISTANCE MATRIX FOR CHECKING
   sprintf(fname,"%s/d-input.dat",outdir);
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<nseq;i++)
       for (j=0;j<nseq;j++)
           fprintf(fp,"%f\n",dist_matrix_pred[i*nseq + j]);
     fclose(fp);

   sprintf(fname,"%s/d-model.dat",outdir);
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<nseq;i++)
       for (j=0;j<nseq;j++)
           fprintf(fp,"%f\n",dist_matrix[i*nseq + j]);
     fclose(fp);

// WRITE LAST DISTANCE MATRIX
   sprintf(fname,"%s/d3-last.dat",outdir);
   //strcpy(fname,"d.dat");
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<nseq;i++)
       for (j=0;j<nseq;j++)
           fprintf(fp,"%f\n",dist_matrix[ind_max[j]*nseq + ind_max[i]]);
     fclose(fp);

// WRITE SCORES
   sprintf(fname,"%s/scores.dat",outdir);
   if ( (fp = fopen(fname,"w")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<NSTEPS;i++)
           fprintf(fp,"%f\n",s[i]);
     fclose(fp);

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Runtime:%f", time_spent);
}
