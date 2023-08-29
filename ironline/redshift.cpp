#include "def.h"
/*
double radius;
double ktt;
double ktkp;
double kyy;
double gg[1];
double ldr[1];
 */
void redshift(double spin, double spin2, double epsilon_r, double epsilon_t, double radius, double ktt, double ktkp, double kyy, double& gg, double& ldr)
{
	double r, r2, r3;
	double dr, temp;
	double cc;
	double def_r,def_t;
	double g00[3], g03[3], g33[3];
	double g001, g031, g331;
	double Omega;
	double uet;
	double uephi;
	double mem;
	double H;
	
	int i;
	
	
	dr   = 0.001*radius;
	temp = radius + dr;
	dr   = temp - radius;
	
	cc = ktkp;
	
	for (i = 0; i <= 2; i++) {
		
		r  = radius + dr*(i - 1);
		r2 = r*r;
		r3 = r2*r;
		
		def_t = epsilon_t/r3;
		def_r = epsilon_r/r3;
		H = (1+def_r)*(1+def_t);
		H = sqrt(H);

		g00[i] = - (1 + def_t)*(1 - 2/r);
		g03[i] = - spin*(H-(1 + def_t)*(1 - 2/r));
		g33[i] = r2 + spin2*(2*H- (1 + def_t)*(1 - 2/r));
		
	}
	
	g001 = 0.5*(g00[2] - g00[0])/dr;
	g031 = 0.5*(g03[2] - g03[0])/dr;
	g331 = 0.5*(g33[2] - g33[0])/dr;
	
	Omega  = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
	
	uet = sqrt(-g00[1] - 2*g03[1]*Omega - g33[1]*Omega*Omega);
	
	gg = uet/(1 - cc*Omega);
	
	uet = 1/uet;
	
	uephi = Omega*uet;
	
	mem = kyy*radius;
	mem = mem/(ktt*uet + ktt*cc*uephi);
	
	if (mem < 0) mem = 0;
	if (mem > 1) mem = 1;
	
	ldr = 0.5 + 0.75*mem;

	return ;
}
