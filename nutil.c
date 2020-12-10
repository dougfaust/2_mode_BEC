// nutil implementation file

#include "math.h"
#include <stdio.h>
#include "nrutil.h"
#include <stdlib.h>
#define NRANSI
#define MAXSTP 10000
#define TINY 1.0e-30
#define NR_END 1
#define FREE_ARG char*

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

double sech( double xx)
{
	double yy,zz;
	if( xx >= 30.0 ) {yy = 0.0; return yy;}
	if( xx <= -30.0) {yy = 0.0; return yy;}
	zz = exp(xx);
	yy = 2.0/(zz + 1.0/zz);
	return yy;
}

// add some Num Rec utility routines

void nrerror(char error_text[])
// std error handler
{	fprintf(stderr, "Run time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...exit to sysyem...\n");
	exit(1);
}


unsigned long *lvector(long nl, long nh)
// allocate an unsigned long vector with range v[nl...nh]
{ unsigned long  *v;
v = (unsigned long *)malloc( (size_t) ((nh-nl + 1 + NR_END)*sizeof(long)));
if(!v) nrerror("allocation failure in vector()");
return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
//free a double vector, orig allocated with lvector()
{free ((FREE_ARG) (v+nl-NR_END));
}

int *ivector(long nl, long nh)
// allocate an unsigned long vector with range v[nl...nh]
{ int  *v;
v = (int *)malloc( (size_t) ((nh-nl + 1 + NR_END)*sizeof(int)));
if(!v) nrerror("allocation failure in vector()");
return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
//free a double vector, orig allocated with lvector()
{free ((FREE_ARG) (v+nl-NR_END));
}

double *vector(long nl, long nh)
// allocate a double vector with range v[nl...nh]
{ double *v;
v = (double *)malloc( (size_t) ((nh-nl + 1 + NR_END)*sizeof(double)));
if(!v) nrerror("allocation failure in vector()");
return v-nl+NR_END;
}

void free_vector(double *v, long nl, long nh)
//free a double vector, orig allocated with vector()
{free ((FREE_ARG) (v+nl-NR_END));
}

double **matrix(long nrl, long nrh, long ncl, long nch)
// allocate  double matrix with subscript range [nrl...nrh][ncl...nch]
{ long i, nrow = nrh-nrl +1, ncol = nch-ncl + 1;
	double **m;

// allocate pointers to rows

	m = (double **)malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if(!m) nrerror("allocation failure 1 in matrix");
	m += NR_END;
	m -= nrl;

// allocate rows and set pointers to them

	m[nrl] = (double *)malloc((size_t)((nrow*ncol + NR_END)*sizeof(double)));
	if(!m[nrl]) nrerror("allocation failure 2 in matrix");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1; i<=nrh;i++) m[i]=m[i-1] + ncol;
// return pointer to array of pointers to rows
	return m;
}

void free_matrix( double **m, long nrl, long nrh, long ncl, long nch)
// frees up storage for double matrix
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m + nrl - NR_END));
}

//		int     decimal,   sign;
//		char    *buffer;
//		int     precision = 10;
//		int     source1 = index;
//		int     source2 = xind;
//		buffer = _ecvt( source1, precision, &decimal, &sign );
//		dc.TextOut(200,100,buffer);
//		buffer = _ecvt( source2, precision, &decimal, &sign );
//		dc.TextOut(200,150,buffer);

