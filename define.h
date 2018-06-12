/*  ȡTi=2Te,������������ΪC3+,������һ�ֲ�����PIC����ó� */
#include <math.h>
/* �������ӵĽṹ */
struct particle 
{
	double z;                                            //λ��
	double vx;                                           //�ٶ�
	double vy; 
	double vz;                                           
	double ener;										//��λeV
	double weig;                                       //Ȩ��
};

/* ����һЩ���� */
static const double pi=3.141592653589793;
static const double k0=1.380658e-23;				//������������
static const double qe=-1.6021865314e-19;			//���ӵ��
static const double qi=1.6021865314e-19;			//���ӵ��
static const double q_C3=4.8065595942e-19;			//C3+���
static const double me=9.1093897e-31;				//��������
static const double m_H=1.674e-27;					//������
static const double m_D=3.348e-27;					//�����
static const double m_C=2.0088e-26;					//̼����
static const double m_W=3.08016e-25;				//������
static const double am_H=1.008;						//H��ԭ������
static const double am_D=2.016;						//D��ԭ������
static const double am_W=183.85;					//�ٵ�ԭ������
static const double am_C=12.01;						//C��ԭ������
static const double z_H=1.0;						//H�ĺ˵����
static const double z_D=1.0;						//D�ĺ˵����
static const double z_W=74.0;						//W�ĺ˵����
static const double z_C=6.0;						//C�ĺ˵����
static const double eps0=8.854187817e-12;			//
static const double Re=0.3;							//���¶�������¶ȱ�ֵ
static const double Te=20.0;						//70����
static const double Ti=Te;
static const double Tc=Ti;							//C3+���¶�
static const double ve=sqrt(-2.0*qe*Te/me);			//sqrt(2.0*qe*Te/me);
static const double Ne=1e19;						//3e18����
static const double B0=5.3;							//2.25����
static const double alpha=(87.0*pi)/180.0;			//�ų���ƫ�����а巨�߷���н�
static const double beta=(60.0*pi)/180.0;			//��λ��
static double bx=sin(alpha)*cos(beta);
static double by=sin(alpha)*sin(beta);
static double bz=-cos(alpha);
static const double lz_plas=0.006;					//�����������򳤶�
static const double lz_sour=0.006;					//Դ���򳤶�
static const double lx=1.0;
static const double ly=1.0;
static const double area=lx*ly;

static const int Nz_plas=1001;							//z����������
static const int Nz_sour=1000;
static const double dz_plas=lz_plas/double(Nz_plas-1);		//z����ռ䲽��
static const double dz_sour=dz_plas;
static const int numb_e=1.0e6;								//ģ����ӵĳ�������
static const int numb_D=1.0e6;								//ģ��D�ĳ�������
static const int numb_C3=0;									//ģ���������ӵĳ�������
static const double weig_e=(lz_plas*Ne)/numb_e;
static const double weig_D=(lz_plas*Ne)/numb_e;
static const double weig_C3=(lz_plas*Ne)/(numb_e*100.0);   //����ģ�����ӵ�Ȩ��
static const int writ_step=500;
static double const U0=0.0;								//-3.0*Te;
static double const vb_D=sqrt((Te+Ti)*qi/m_D);			//뮵���������
static double const vb_C=sqrt((Te+Tc)*qi/m_C);			//C���ӵ���������

static const double dt=1.0e-12;						//ʱ�䲽��
static const double freq_D=qe*B0/m_D;				//����ӻ�����Ƶ��
static const double radi_D=sqrt(3.0*m_D*Ti/qe);		//����ӻ����뾶
static const double deby_leng=sqrt(-eps0*Te*Ti/((Te+Ti)*Ne*qe));

/* �ٰ�ԭ�ӵı���������(eV)���ܶ�(g/cm3) */
static const double ener_sur=8.7;
static const double dens_targ=19.35;
static const double NA=6.0221367e23;                  //����٤���޳���

void diag(void);
void init(particle *ptr_e,particle *ptr_D,particle *ptr_C3);																	//��ʼ��
void maxw(double *vx,double *vy,double *vz,double kt,double mass);																//Maxwell�ٶȷֲ�
void pic(particle *ptr,int numb,double q,double *dens);																			//�����Ʒ��ҷ�
void grid_pic(particle *prt_e,particle *ptr_D,particle *ptr_C3,int *numb,double (*dens)[Nz_plas],double *char_dens);			
void pois(double *char_dens,int n,double * pote);																				//��Ⲵ�ɷ���
void solvfiel(double *pote,int n,double *fiel);																					//���糡
void setv(particle *ptr,int numb,double *fiel,double q,double m);																//��⣨-0.5dt)ʱ���ٶ�
void move(particle *ptr,int *numb,double *fiel,double q,double m,double *flux,double *ener_flux);									//�ƶ�����