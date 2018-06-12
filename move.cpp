#include <stdio.h>
#include <math.h>
#include "define.h"
void move(particle *ptr,int *numb,double *fiel,double q,double m,double *flux,double *ener_flux)
{
	int i,j,ip;
	double tt[3],u[3],tt1,u1[3],u2[3],u3[3],tan_b;									//���Boris�㷨�еĸ��ٶȷ�����Ϣ
	double z,ez,az;
	double l1,l2;
	double vx,vy,vz;
	double weig,ener;
	double cos_thet;																//���䷽����а巨�߷���н�cosֵ

	j=0;																			//ͳ�ƻ��ڼ��������ڵ�������
	for(i=0;i<2;i++)
	{
		flux[i]=0.0;
		ener_flux[i]=0.0;
	}

	for(i=0;i<*numb;i++)
	{
		z=ptr[i].z;
		vx=ptr[i].vx;
		vy=ptr[i].vy;
		vz=ptr[i].vz;
		weig=ptr[i].weig;
		/* ���糡 */
		ip=(int)(z/dz_plas);
		l2=z/dz_plas-ip;
		l1=1.0-l2;
		ez=fiel[ip]*l1+fiel[ip+1]*l2;
		az=q*ez/m;
		/* Boris �㷨 J. P. Verboncoeur, Plasma Phys Contr F 47, A231 (2005)  ���ڴų������ 
		tan_b=tan(q*B0*dt/(2.0*m));
		tt[0]=tan_b*bx;
		tt[1]=tan_b*by;
		tt[2]=tan_b*bz;
		tt1=2.0/(1.0+tt[0]*tt[0]+tt[1]*tt[1]+tt[2]*tt[2]);
		u[0]=vx;
		u[1]=vy;
		u[2]=vz;
		u1[0]=u[0];
		u1[1]=u[1];
		u1[2]=u[2]+0.5*az*dt;
		u2[0]=u1[0]+(u1[1]*tt[2]-u1[2]*tt[1]);
		u2[1]=u1[1]+(u1[2]*tt[0]-u1[0]*tt[2]);
		u2[2]=u1[2]+(u1[0]*tt[1]-u1[1]*tt[0]);
		u3[0]=u1[0]+tt1*(u2[1]*tt[2]-u2[2]*tt[1]);
		u3[1]=u1[1]+tt1*(u2[2]*tt[0]-u2[0]*tt[2]);
		u3[2]=u1[2]+tt1*(u2[0]*tt[1]-u2[1]*tt[0]);
		u[0]=u3[0];
		u[1]=u3[1];
		u[2]=u3[2]+0.5*az*dt;
		vx=u[0];
		vy=u[1];
		vz=u[2];   */
		vz=vz+az*dt;

		z=z+vz*dt;

		/*�ж��Ƿ񵽴�߽�,���������ձ߽����� */
		if((z>0.0)&&(z<lz_plas))
		{
			ptr[j].z=z;
			ptr[j].vx=vx;
			ptr[j].vy=vy;
			ptr[j].vz=vz;
			ener=0.5*m*(vx*vx+vy*vy+vz*vz)/qi;
			ptr[j].ener=ener;
			ptr[j].weig=weig;
			j++;
		}
		else if(z<=0.0)
		{
			ener=0.5*m*(vx*vx+vy*vy+vz*vz)/qi;												//��λeV,����Ҫ����ת��
			cos_thet=fabs(vz)/sqrt(vx*vx+vy*vy+vz*vz);
			flux[0]=(flux[0]+weig)*cos_thet;
			ener_flux[0]=(ener_flux[0]+weig*ener)*cos_thet;

		}
		else if(z>=lz_plas)
		{
			ener=0.5*m*(vx*vx+vy*vy+vz*vz)/qi;												//��λeV,����Ҫ����ת��
			cos_thet=fabs(vz)/sqrt(vx*vx+vy*vy+vz*vz);
			flux[1]=(flux[1]+weig)*cos_thet;
			ener_flux[1]=(ener_flux[1]+weig*ener)*cos_thet;

		}
	}

	*numb=j;
	for(i=0;i<2;i++)
	{
		flux[i]=flux[i]/(area*dt);
		ener_flux[i]=ener_flux[i]*qi/(area*dt);
	}
}
