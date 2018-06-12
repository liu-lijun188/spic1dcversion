/*  取Ti=2Te,入射杂质离子为C3+,能量均一分布，由PIC计算得出 */
#include <math.h>
/* 定义离子的结构 */
struct particle 
{
	double z;                                            //位置
	double vx;                                           //速度
	double vy; 
	double vz;                                           
	double ener;										//单位eV
	double weig;                                       //权重
};

/* 定义一些常量 */
static const double pi=3.141592653589793;
static const double k0=1.380658e-23;				//玻耳兹曼常数
static const double qe=-1.6021865314e-19;			//电子电荷
static const double qi=1.6021865314e-19;			//离子电荷
static const double q_C3=4.8065595942e-19;			//C3+电荷
static const double me=9.1093897e-31;				//电子质量
static const double m_H=1.674e-27;					//氢质量
static const double m_D=3.348e-27;					//氘质量
static const double m_C=2.0088e-26;					//碳质量
static const double m_W=3.08016e-25;				//钨质量
static const double am_H=1.008;						//H的原子质量
static const double am_D=2.016;						//D的原子质量
static const double am_W=183.85;					//钨的原子质量
static const double am_C=12.01;						//C的原子质量
static const double z_H=1.0;						//H的核电荷数
static const double z_D=1.0;						//D的核电荷数
static const double z_W=74.0;						//W的核电荷数
static const double z_C=6.0;						//C的核电荷数
static const double eps0=8.854187817e-12;			//
static const double Re=0.3;							//钨温度与电子温度比值
static const double Te=20.0;						//70测试
static const double Ti=Te;
static const double Tc=Ti;							//C3+的温度
static const double ve=sqrt(-2.0*qe*Te/me);			//sqrt(2.0*qe*Te/me);
static const double Ne=1e19;						//3e18测试
static const double B0=5.3;							//2.25测试
static const double alpha=(87.0*pi)/180.0;			//磁场与偏滤器靶板法线方向夹角
static const double beta=(60.0*pi)/180.0;			//方位角
static double bx=sin(alpha)*cos(beta);
static double by=sin(alpha)*sin(beta);
static double bz=-cos(alpha);
static const double lz_plas=0.006;					//等离子体区域长度
static const double lz_sour=0.006;					//源区域长度
static const double lx=1.0;
static const double ly=1.0;
static const double area=lx*ly;

static const int Nz_plas=1001;							//z方向网格数
static const int Nz_sour=1000;
static const double dz_plas=lz_plas/double(Nz_plas-1);		//z方向空间步长
static const double dz_sour=dz_plas;
static const int numb_e=1.0e6;								//模拟电子的超粒子数
static const int numb_D=1.0e6;								//模拟D的超粒子数
static const int numb_C3=0;									//模拟杂质粒子的超粒子数
static const double weig_e=(lz_plas*Ne)/numb_e;
static const double weig_D=(lz_plas*Ne)/numb_e;
static const double weig_C3=(lz_plas*Ne)/(numb_e*100.0);   //各种模拟粒子的权重
static const int writ_step=500;
static double const U0=0.0;								//-3.0*Te;
static double const vb_D=sqrt((Te+Ti)*qi/m_D);			//氘的离子声速
static double const vb_C=sqrt((Te+Tc)*qi/m_C);			//C离子的离子声速

static const double dt=1.0e-12;						//时间步长
static const double freq_D=qe*B0/m_D;				//氘离子回旋角频率
static const double radi_D=sqrt(3.0*m_D*Ti/qe);		//氘离子回旋半径
static const double deby_leng=sqrt(-eps0*Te*Ti/((Te+Ti)*Ne*qe));

/* 钨靶原子的表面束缚能(eV)及密度(g/cm3) */
static const double ener_sur=8.7;
static const double dens_targ=19.35;
static const double NA=6.0221367e23;                  //阿伏伽德罗常数

void diag(void);
void init(particle *ptr_e,particle *ptr_D,particle *ptr_C3);																	//初始化
void maxw(double *vx,double *vy,double *vz,double kt,double mass);																//Maxwell速度分布
void pic(particle *ptr,int numb,double q,double *dens);																			//粒子云分室法
void grid_pic(particle *prt_e,particle *ptr_D,particle *ptr_C3,int *numb,double (*dens)[Nz_plas],double *char_dens);			
void pois(double *char_dens,int n,double * pote);																				//求解泊松方程
void solvfiel(double *pote,int n,double *fiel);																					//求解电场
void setv(particle *ptr,int numb,double *fiel,double q,double m);																//求解（-0.5dt)时刻速度
void move(particle *ptr,int *numb,double *fiel,double q,double m,double *flux,double *ener_flux);									//推动粒子