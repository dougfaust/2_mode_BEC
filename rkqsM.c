//#include "stdafx.h"
//#include "BoxBEC.h"
//#include "BoxBECDlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

#include <math.h>
#include "nrutil.h"

#define NRANSI
#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void rkqsM(double y[],double dydx[],int n,double *x,double htry, double eps, double yscal[], double *hdid, double *hnext, void(*derivs2)(double [], double [], int))
{
		void rkckM(double y[], double dydx[], int n, double x, double h, double yout[], double yerr[], void (*derivs2)(double [], double [], int));
		
	int i; int num=1;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=vector(0,(long)(n-1)/2);
	ytemp=vector(0,(long)n);
	
	h=htry;
	for(;;) 
	{
		rkckM(y, dydx, n, *x, h, ytemp, yerr, derivs2);
		errmax=0.0;
		for (i=0;i<=(n-1)/2;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
		xnew=(*x)+h;
		if (xnew == *x) nrerror("stepsize underflow in rkqs");
	}
	if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
	else *hnext=5.0*h;
	*x += (*hdid=h);
	for (i=0;i<=n;i++) {y[i]=ytemp[i];}
	free_vector(ytemp,0,(long)n);
	free_vector(yerr,0,(long)(n-1)/2);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI






