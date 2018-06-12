#include <stdio.h>
#include "define.h"
// n=Nz_plas
void solvfiel(double *pote,int n,double *fiel)
{
	int i;
	fiel[0]=-(-3.0*pote[0]+4.0*pote[1]-1.0*pote[2])/(2.0*dz_plas);
	for(i=1;i<n-1;i++)
		fiel[i]=-(pote[i+1]-pote[i-1])/(2.0*dz_plas);
	fiel[n-1]=-(1.0*pote[n-3]-4.0*pote[n-2]+3.0*pote[n-1])/(2.0*dz_plas);

}