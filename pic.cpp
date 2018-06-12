#include <stdio.h>
#include <math.h>
#include "define.h"
/*  *ptr_e,ptr_D,ptr_C3分别代表指向三种类型的指针
	* *dens[]代表三种粒子的密度
	*char_dens代表电荷密度		*/
void grid_pic(particle *prt_e,particle *ptr_D,particle *ptr_C3,int *numb,double (*dens)[Nz_plas],double *char_dens)
{
	int i;

	pic(prt_e,numb[0],qe,dens[0]);
	pic(ptr_D,numb[1],qi,dens[1]);
	pic(ptr_C3,numb[2],q_C3,dens[2]);
	for(i=0;i<Nz_plas;i++)
		char_dens[i]=dens[0][i]+dens[1][i]+dens[2][i];
	for(i=0;i<Nz_plas;i++)
	{
		dens[0][i]/=qe;
		dens[1][i]/=qi;
		dens[2][i]/=q_C3;
	}
}

/* 云分室法求解电荷密度 */
void pic(particle *ptr,int numb,double q,double *dens)
{
	int i,ip;
	double l1,l2;
	double q_real,z;

	for(i=0;i<Nz_plas;i++)
		dens[i]=0.0;

	for(i=0;i<numb;i++)
	{
		q_real=ptr[i].weig*q;
		z=ptr[i].z;
		ip=(int)(z/dz_plas);
		l2=z/dz_plas-ip;
		l1=1.0-l2;
		dens[ip]+=q_real*l1;
		dens[ip+1]+=q_real*l2;
	}

	for(i=1;i<Nz_plas-1;i++)
		dens[i]/=dz_plas;
	dens[0]/=(0.5*dz_plas);				
	dens[Nz_plas-1]/=(0.5*dz_plas);
}