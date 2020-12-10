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

void rkqsW(double y[],double dydx[],int n,double *x,double htry
		  ,double eps, double yscal[],double *hdid, double *hnext,
		  double ke_x[], double scale_factor, unsigned long nvar_x,
		  void(*derivs)(double, double [],double [],double[],double,unsigned long))
{
		void rkckW(double y[], double dydx[], int n,
		  double x, double h, double yout[],
	double yerr[], double ke_x[], double scale_factor, unsigned long nvar_x,
	void (*derivs)( double, double [], double [],double[],double,unsigned long) );
		
	int i; int num=1;
	double errmax,h,htemp,xnew,*yerr,*ytemp;

	yerr=vector(1,(long)n/2);
	ytemp=vector(1,(long)n);
	
	h=htry;
	for(;;) 
	{
		rkckW(y,dydx,n,*x,h,ytemp,yerr,ke_x, scale_factor, nvar_x,derivs);
		errmax=0.0;
		for (i=1;i<=n/2;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
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
	for (i=1;i<=n;i++) {y[i]=ytemp[i];}
	free_vector(ytemp,1,(long)n);
	free_vector(yerr,1,(long)n/2);
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI






