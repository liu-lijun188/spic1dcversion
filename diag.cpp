#include "define.h"
#include <stdio.h>
#include <math.h>
void diag(void)
{
	double vD;
	double lamd_e,lamd_D;
	double w_pe,peri_e;
	double wpe_dt,vdt_dz;
	vD=sqrt(2.0*Ti*qi/m_D);				//氘离子热速度
	lamd_e=me*ve/fabs(qe*B0);			//电子回旋半径
	lamd_D=m_D*vD/fabs(qi*B0);			//氘离子回旋半径
	w_pe=sqrt(Ne*qe*qe/(eps0*me));		//电子震荡频率
	peri_e=2.0*pi/w_pe;
	wpe_dt=w_pe*dt;
	vdt_dz=ve*dt/dz_plas;
	printf("ne=%10.6e,Te=%10.6f\n",Ne,Te);
	printf("ve=%10.6e,v_D=%10.6e\n",ve,vD);
	printf("德拜长度=%10.6e,空间步长=%10.6e\n",deby_leng,dz_plas);
	printf("电子振荡周期=%10.6e,时间步长=%10.6e\n",peri_e,dt);
	printf("wpe*dt=%10.6f,v*dt/dz=%10.6f\n",wpe_dt,vdt_dz);
	
}