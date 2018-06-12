#include "define.h"
#include <stdio.h>
#include <math.h>
void diag(void)
{
	double vD;
	double lamd_e,lamd_D;
	double w_pe,peri_e;
	double wpe_dt,vdt_dz;
	vD=sqrt(2.0*Ti*qi/m_D);				//��������ٶ�
	lamd_e=me*ve/fabs(qe*B0);			//���ӻ����뾶
	lamd_D=m_D*vD/fabs(qi*B0);			//����ӻ����뾶
	w_pe=sqrt(Ne*qe*qe/(eps0*me));		//������Ƶ��
	peri_e=2.0*pi/w_pe;
	wpe_dt=w_pe*dt;
	vdt_dz=ve*dt/dz_plas;
	printf("ne=%10.6e,Te=%10.6f\n",Ne,Te);
	printf("ve=%10.6e,v_D=%10.6e\n",ve,vD);
	printf("�°ݳ���=%10.6e,�ռ䲽��=%10.6e\n",deby_leng,dz_plas);
	printf("����������=%10.6e,ʱ�䲽��=%10.6e\n",peri_e,dt);
	printf("wpe*dt=%10.6f,v*dt/dz=%10.6f\n",wpe_dt,vdt_dz);
	
}