// sinft program
// Calculates the sine transform of a set of n real-valued data points stored in array 
// y[1..n].  The number n must be a power of 2.  On exit y is replaced by its tranform.
// This program also calculates the inverse sine transform. but the output array should
// be multiplied by 2/n.

//#include "stdafx.h"
//#include "BoxBEC.h"
//#include "BoxBECDlg.h"
#include <math.h>

void sinft(double y[], int n)

{
	void realft(double data[], unsigned long n, int isign);
	int j, n2=n+2;
	double sum, y1, y2;
	double theta, wi=0.0, wr=1.0, wpi, wpr, wtemp;
	theta=3.14159265358979/(double)n;
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	y[1]=0.0;
	for (j=2;j<=(n>>1)+1;j++) 	
	{
	    wr=(wtemp=wr)*wpr-wi*wpi+wr;   // Calculate the sine for the auxillary array.
	    wi=wi*wpr+wtemp*wpi+wi;	   // The cosine is needed to scontinue the recurrence.
	    y1=wi*(y[j]+y[n2-j]);	   // Construct the auxiliary array.
	    y2=0.5*(y[j]-y[n2-j]);          
	    y[j]=y1+y2;			   // Term j and N-j are related.
	    y[n2-j]=y1-y2;
	}
	realft(y,n,1);			   // Transform the auxiliary array.
	y[1]*=0.5;			   // Initialize the sum used for odd terms below.
	sum=y[2]=0.0;
	for (j=1;j<=n-1;j+=2) 	{
	    sum += y[j];		   // Even terms determined directly. 
	   y[j]=y[j+1];
		y[j+1]=sum;			   // Odd terms determined by this running sum.
	}
}


