#include <stdio.h>
#include <stdlib.h>
#include "define.h"
//n=Nz_plas
/* 求解泊松方程 第一类边界条件 */
void pois(double *char_dens,int n,double * pote)
{
	double *a,*b,*c,*f,*e,*d;
	int i;
	a=(double *)malloc(n*sizeof(double));
	b=(double *)malloc(n*sizeof(double));
	c=(double *)malloc(n*sizeof(double));
	f=(double *)malloc(n*sizeof(double));
	e=(double *)malloc(n*sizeof(double));
	d=(double *)malloc(n*sizeof(double));

	pote[0]=U0;
	pote[n-1]=U0;
	f[1]=-char_dens[1]*dz_plas*dz_plas/eps0-pote[0];
	f[n-2]=-char_dens[n-2]*dz_plas*dz_plas/eps0-pote[n-1];
	for(i=2;i<n-2;i++)
		f[i]=-char_dens[i]*dz_plas*dz_plas/eps0;
	for(i=0;i<n;i++)
	{
		a[i]=1.0;
		b[i]=-2.0;
		c[i]=1.0;
	}
	e[1]=c[1]/b[1];
	d[1]=f[1]/b[1];
	for(i=2;i<n-2;i++)
	{
		e[i]=c[i]/(b[i]-a[i]*e[i-1]);
		d[i]=(f[i]-a[i]*d[i-1])/(b[i]-a[i]*e[i-1]);
	}
	pote[n-2]=(f[n-2]-a[n-2]*d[n-3])/(b[n-2]-a[n-2]*e[n-3]);
	for(i=n-3;i>0;i--)
		pote[i]=d[i]-e[i]*pote[i+1];
	free(a);
	free(b);
	free(c);
	free(f);
	free(e);
	free(d);
}