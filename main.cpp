#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "define.h"
#define NUMBMAX 3000000					
#define WRITSPACSTEP 1			//空间每隔多少步写入
#define SIMUSTEP 25000			//模拟总步数
#define SHOWSTEP 500			//时间上每隔多少步求平均值
#define SHOWNUMB 13				
int main(void)
{
	static particle part_e[NUMBMAX],part_D[NUMBMAX],part_C3[NUMBMAX];			//三种模拟粒子
	int numb_edc[3];															//存放三种粒子的粒子数
	double dens[3][Nz_plas];													//存放三种粒子在各格点的密度
	double char_dens[Nz_plas],pote[Nz_plas],fiel[Nz_plas];						//存放各格点的电荷密度、电势、电场
	double flux[3][2],ener_flux[3][2],flux_aver[3][2],ener_flux_aver[3][2];		//统计入射到边界的粒子流及能流
	FILE *fp_init,*fp_aver,*fp_flux_aver,*fp_ener_flux_aver,*fp_dura;
	double dens_aver[3][Nz_plas],char_dens_aver[Nz_plas],pote_aver[Nz_plas],fiel_aver[Nz_plas]; //存放平均值
	int i,j,k;
	FILE *fp_who[SHOWNUMB];
	int show_who[SHOWNUMB]={1000,2000,3000,5000,7000,10000,15000,25000,35000,45000,60000,80000,100000};
	clock_t star,fini;
	int hour,minu;
	star=clock();
	/* 初始化 */
	numb_edc[0]=numb_e;
	numb_edc[1]=numb_D;
	numb_edc[2]=numb_C3;
	for(i=0;i<Nz_plas;i++)
	{
		for(j=0;j<3;j++)
			dens_aver[j][i]=0.0;
		char_dens_aver[i]=0.0;
		pote_aver[i]=0.0;
		fiel_aver[i]=0.0;
	}
	for(i=0;i<3;i++)
	{
		for(j=0;j<2;j++)
		{
			flux_aver[i][j]=0.0;
			ener_flux_aver[i][j]=0.0;
		}
	}
	fp_init=fopen("intial.dat","w");
	fp_aver=fopen("aver.dat","w");
	fp_flux_aver=fopen("flux_aver.dat","w");
	fp_ener_flux_aver=fopen("ener_flux_aver.dat","w");
	fp_dura=fopen("duration.dat","w");
	fp_who[0]=fopen("1ns.dat","w");
	fp_who[1]=fopen("2ns.dat","w");
	fp_who[2]=fopen("3ns.dat","w");
	fp_who[3]=fopen("5ns.dat","w");
	fp_who[4]=fopen("7ns.dat","w");
	fp_who[5]=fopen("10ns.dat","w");
	fp_who[6]=fopen("15ns.dat","w");
	fp_who[7]=fopen("25ns.dat","w");
	fp_who[8]=fopen("35ns.dat","w");
	fp_who[9]=fopen("45ns.dat","w");
	fp_who[10]=fopen("60ns.dat","w");
	fp_who[11]=fopen("80ns.dat","w");
	fp_who[12]=fopen("100ns.dat","w");

	diag();
	srand((unsigned int)time(0));
	printf("Initial......\n");
	init(part_e,part_D,part_C3);												//初始化
	grid_pic(part_e,part_D,part_C3,numb_edc,dens,char_dens);					//PIC
	pois(char_dens,Nz_plas,pote);												//求解电势
	solvfiel(pote,Nz_plas,fiel);												//求解电场
	/*  测试初始化是否正确 */
	for(i=0;i<Nz_plas;i+=WRITSPACSTEP)
		fprintf(fp_init,"%10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",dens[0][i],dens[1][i],dens[2][i],char_dens[i],pote[i],fiel[i]);
	/* 求解-0.5dt时刻速度 */
	setv(part_e,numb_edc[0],fiel,qe,me);
	setv(part_D,numb_edc[1],fiel,qi,m_D);
	setv(part_C3,numb_edc[2],fiel,q_C3,m_C);

	/* 开始循环求解 */
	for(i=1;i<=SIMUSTEP;i++)
	{
		if(i%SHOWSTEP==0)
			printf("i=%d\n",i);
		/* 推动粒子 */
		move(part_e,&numb_edc[0],fiel,qe,me,flux[0],ener_flux[0]);         //电子
		move(part_D,&numb_edc[1],fiel,qi,m_D,flux[1],ener_flux[1]);			//离子
		move(part_C3,&numb_edc[2],fiel,q_C3,m_C,flux[2],ener_flux[2]);		//杂质粒子

		grid_pic(part_e,part_D,part_C3,numb_edc,dens,char_dens);
		pois(char_dens,Nz_plas,pote);
		solvfiel(pote,Nz_plas,fiel);

		/* 统计各物理量的平均值 */
		for(j=0;j<Nz_plas;j++)							
		{
			for(k=0;k<3;k++)
				dens_aver[k][j]+=dens[k][j];
			char_dens_aver[j]+=char_dens[j];
			pote_aver[j]+=pote[j];
			fiel_aver[j]+=fiel[j];
		}
		for(j=0;j<3;j++)
			for(k=0;k<2;k++)
			{
				flux_aver[j][k]+=flux[j][k];
				ener_flux_aver[j][k]+=ener_flux[j][k];
			}

		if(i%SHOWSTEP==0)
		{
			for(j=0;j<Nz_plas;j++)
			{
				for(k=0;k<3;k++)
					dens_aver[k][j]/=(double)SHOWSTEP;
				char_dens_aver[j]/=(double)SHOWSTEP;
				pote_aver[j]/=(double)SHOWSTEP;
				fiel_aver[j]/=(double)SHOWSTEP;
			}
			for(j=0;j<3;j++)
				for(k=0;k<2;k++)
				{
					flux_aver[j][k]/=(double)SHOWSTEP;
					ener_flux_aver[j][k]/=(double)SHOWSTEP;
				}

			//输出平均电势等信息
			fprintf(fp_aver,"time=%10.6e\n",i*dt);
			for(j=0;j<Nz_plas;j+=WRITSPACSTEP)
			{
				fprintf(fp_aver,"%15.6e %15.6e  %15.6e  %15.6e  %15.6e  %15.6e  %15.6e\n",j*dz_plas,dens_aver[0][j],dens_aver[1][j],dens_aver[2][j],char_dens_aver[j],pote_aver[j],fiel_aver[j]);
			}
			//输出达到边界的平均粒子流、能流
			for(j=0;j<3;j++)
				for(k=0;k<2;k++)
				{
					fprintf(fp_flux_aver,"%15.6e ",flux_aver[j][k]);
					fprintf(fp_ener_flux_aver,"%15.6e ",ener_flux_aver[j][k]);
				}
			fprintf(fp_flux_aver,"\n");
			fprintf(fp_ener_flux_aver,"\n");

			//输出瞬时电势等信息
			for(j=0;j<SHOWNUMB;j++)
				if(i==show_who[j])
					for(k=0;k<Nz_plas;k++)
						fprintf(fp_who[j],"%12.4e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n",k*dz_plas,dens[0][k],dens[1][k],dens[2][k],char_dens[k],pote[k],fiel[k]);

			//清零
			for(j=0;j<Nz_plas;j++)
			{
				for(k=0;k<3;k++)
					dens_aver[k][j]=0.0;
				char_dens_aver[j]=0.0;
				pote_aver[j]=0.0;
				fiel_aver[j]=0.0;
			}

			for(j=0;j<3;j++)
				for(k=0;k<2;k++)
				{
					flux_aver[j][k]=0.0;
					ener_flux_aver[j][k]=0.0;
				}
		}

	}

	fclose(fp_init);
	fclose(fp_aver);
	fclose(fp_flux_aver);
	fclose(fp_ener_flux_aver);
	for(i=0;i<13;i++)
		fclose(fp_who[i]);
	fini=clock();
	hour=((fini-star)/CLOCKS_PER_SEC)/3600;
	minu=((fini-star)/CLOCKS_PER_SEC)/60-hour*60;
	fprintf(fp_dura,"This program runs  %10d hours %10d minutes\n",hour,minu);
	fclose(fp_dura);
	return 0;
}