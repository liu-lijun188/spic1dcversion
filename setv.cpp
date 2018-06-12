#include "define.h"
#include <stdio.h>
void setv(particle *ptr,int numb,double *fiel,double q,double m)
{
	double vz,az,z,ez;
	int ip,i;
	double l1,l2;

	for(i=0;i<numb;i++)
	{
		z=ptr[i].z;
		vz=ptr[i].vz;
		ip=(int)(z/dz_plas);
		l2=z/dz_plas-ip;
		l1=1.0-l2;
		ez=fiel[ip]*l1+fiel[ip+1]*l2;
		az=q*ez/m;
		vz-=(0.5*az*dt);
		ptr[i].vz=vz;
	}

}