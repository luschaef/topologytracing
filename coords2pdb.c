
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include <unistd.h>
#define MAX_STRING_SIZE 9999

int arg_count;
int j,i;
int n;
double *p;
FILE *fp;
char fname[MAX_STRING_SIZE];

char coordfile[MAX_STRING_SIZE];
char pdbfile[MAX_STRING_SIZE];

int main (int argc, char *argv[])
{
//=============================
// PARSE COMMAND LINE ARGUMENTS
// ============================

  arg_count = 1;
  while ( arg_count < argc ){

   
     if ( ! strcmp(argv[arg_count],"-i") )
     {
         strcpy(coordfile,argv[arg_count+1]);
         arg_count+=2;

     } else if ( ! strcmp(argv[arg_count],"-o") )
     {
         strcpy(pdbfile,argv[arg_count+1]);
         arg_count+=2;
	
     } else if ( ! strcmp(argv[arg_count], "-n") )
     {
	 n = atoi(argv[arg_count+1]);
         arg_count+=2;
     }

  }
// ===================
// READ COORDINATES
// ===================
     if ( (p= (double *) calloc(3*n,sizeof(double)) ) ==  NULL ) {
             fprintf(stderr,"ERROR: could not allocate memory \n");
             exit(1);
     }
     strcpy(fname,coordfile);
     if ( (fp = fopen(fname,"r")) == NULL)
     {
       fprintf(stderr,"ERROR: Can't open file : %s\n",fname);
       exit(1);
     }
     for (i=0;i<n;i++)
           fscanf(fp,"%lf %lf %lf",&p[i*3+0],&p[i*3+1],&p[i*3+2]);

     fclose(fp);

// =========
// WRITE PDB 
// =========
char    atom_name[]="CA";
char    residue_name[]="ALA";
int     position_labels=1;
double  temp_factor=1.00;
char    chain_id[]="A";
double  x1,y1,z1;


fp = fopen(pdbfile,"w");

for (j=0;j<n;j++) {
   x1= p[j*3+0];
   y1= p[j*3+1];
   z1= p[j*3+2];
   fprintf(fp,"ATOM  %5i  %-4s%-4s%5i    %8.3f%8.3f%8.3f %5.2f %5.2f   %4s\n",j, \
            atom_name,residue_name,position_labels++,x1,y1,z1,  \
            temp_factor,temp_factor,chain_id);
}
fprintf(fp,"END \n");
fclose(fp);

}
