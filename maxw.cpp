#include <stdlib.h>
#include <math.h>
#include "define.h"
/* James Joseph Szabo.Fully Kinetic Numerical Modeling of a Plasma Thruster.2011,p149 */

void maxw(double *vx,double *vy,double *vz,double kt,double mass)
{
	double ran1,vt,v;
	double theta,psi;

	vt=sqrt(2.0*kt*qi/mass);
	ran1=rand()/(double)RAND_MAX;
	v=vt*sqrt(-log(ran1));
	ran1=rand()/(double)RAND_MAX;
	theta=2.0*pi*ran1;
	ran1=rand()/(double)RAND_MAX;
	psi=2.0*pi*ran1;

	*vx=v*sin(theta);
	*vy=v*cos(theta);
	*vz=v*sin(psi);

}