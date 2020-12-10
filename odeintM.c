

#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define MAXSTP 5000
#define TINY 1.0e-30
#define NR_END 1
#define FREE_ARG char*

void odeintM(double ystart[], int nvarr, double x1, double x2, double eps, double h1, double hmin, int *nok, int *nbad,
	     void (*derivs2)(double [], double [], int),
	void (*rkqsM)(double [], double [], int, double *, double, double, double [], double *, double *, void (*)(double [], double [], int)))
{
	 
  /* debug 
  FILE *fp2;
  char *f2="derivs.out";
  fp2=fopen(f2, "w");
  //debug  */	 
	 
	int nstp,i;
	double x,hnext,hdid,h;
	double *yscal,*y,*dydx;
	
	yscal=vector(0,(long)nvarr);
	y=vector(0,(long)nvarr);
	dydx=vector(0,(long)nvarr);
	
	x=x1;
	h=SIGN(h1,x2-x1);
	*nok = (*nbad) =  0;
	
	for (i=0;i<=nvarr;i++) {y[i]=ystart[i]; }
	
	for (nstp=1;nstp<=MAXSTP;nstp++) {
		
		(*derivs2)(y,dydx,nvarr);

    /* debug 
    for(i=1;i<nvarr;i++){
      fprintf(fp2, "%d  %d  %20.12f\n", nstp, i, dydx[i]);
    } 
    //debug */


		for (i=0;i<=nvarr;i++)
		{yscal[i]=sqrt(y[i]*y[i])+sqrt(dydx[i]*dydx[i]*h)+TINY;
		if(yscal[i]<=0.001) yscal[i] = 0.001;}
				
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		
		(*rkqsM)(y, dydx, nvarr, &x, h, eps, yscal, &hdid, &hnext, derivs2);
		// printf("in odeintM, hdid and hnext = %20.8e,  %20.8e \n", hdid, hnext);
	
		if (hdid == h) ++(*nok); else ++(*nbad);
		if ((x-x2)*(x2-x1) >= 0.0) 
		{
			for (i=0;i<=nvarr;i++) 
				{ystart[i]=y[i];}
			
			free_vector(dydx,0,(long)nvarr);
			free_vector(y,0,(long)nvarr);
			free_vector(yscal,0,(long)nvarr);
			
			return;
		}
		if (fabs(hnext) <= hmin) nrerror("Step size too small in odeint");
		h=hnext;
	}
  /* debug 
  fclose(fp2);
  //debug */

	nrerror("Too many steps in routine odeint");
}
#undef MAXSTP
#undef TINY
#undef NRANSI

