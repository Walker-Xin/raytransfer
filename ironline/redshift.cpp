#include "def.h"
/*
double radius;
double ktt;
double ktkp;
double kyy;
double gg[1];
double ldr[1];
 */
void redshift_alt(double spin, double spin2, double epsilon_r, double epsilon_t, double radius, double ktt, double ktkp, double kyy, double& gg, double& ldr)
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

void redshift(double spin, double defpar, double r, double th, double ktkp, double &gg)
{
double t8 = pow(spin,2);
double t7 = pow(r,2);
double t9 = t7 + t8;
double t21 = pow(r,3);
double t12 = sin(th);
double t13 = pow(t12,2);
double t22 = defpar + t21;
double t4 = pow(r,-3);
double t5 = defpar*t4;
double t6 = 1 + t5;
double t11 = t6*t9;
double t14 = -(t13*t8);
double t15 = t11 + t14;
double t26 = defpar*t7;
double t27 = pow(r,5);
double t28 = -2*t27;
double t29 = 3*defpar*t8;
double t30 = t26 + t28 + t29;
double t17 = cos(th);
double t18 = pow(t17,2);
double t19 = t18*t8;
double t20 = t19 + t7;
double t40 = pow(t9,2);
double t42 = -2 + r;
double t43 = r*t42;
double t44 = t43 + t8;
double t23 = t22*t9;
double t24 = -(t13*t21*t8);
double t25 = t23 + t24;
double t57 = pow(r,4);
double t58 = 2*t57;
double t59 = defpar*t8;
double t60 = t26 + t58 + t59;
double t61 = 2*defpar;
double t62 = 5*t21;
double t63 = 3*r*t8;
double t64 = -3*r*t13*t8;
double t65 = t61 + t62 + t63 + t64;
double t66 = 2*t20*t60*t65*t7;
double t67 = -2*t25*t60*t7;
double t68 = 2*defpar*r;
double t69 = 8*t21;
double t70 = t68 + t69;
double t71 = -(r*t20*t25*t70);
double t72 = -3*t20*t25*t60;
double t73 = t66 + t67 + t71 + t72;
double t32 = -r;
double t33 = 1 + t32;
double t31 = -(t22*t30*t9);
double t34 = pow(r,7);
double t35 = t13*t33*t34*t8;
double t36 = t31 + t35;
double t37 = t20*t25*t36;
double t38 = pow(r,11);
double t39 = pow(t6,2);
double t41 = t39*t40;
double t45 = -(t13*t44*t8);
double t46 = t41 + t45;
double t47 = t15*t38*t46;
double t48 = pow(t22,2);
double t49 = t40*t48;
double t50 = pow(r,6);
double t51 = -(t13*t44*t50*t8);
double t52 = t49 + t51;
double t53 = t20*t30*t52;
double t54 = t37 + t47 + t53;
double t88 = t14 + t43 + t8;
double t16 = pow(t15,3);
double t55 = 1/t54;
double t56 = pow(t25,-3);
double t74 = -(spin*t13*t56*t7*t73);
double t75 = pow(t12,4);
double t76 = pow(t25,-6);
double t77 = pow(t73,2);
double t78 = t57*t75*t76*t77*t8;
double t79 = pow(r,-14);
double t80 = pow(t15,-6);
double t81 = 2 + t32;
double t82 = r*t81;
double t83 = -t8;
double t84 = t13*t8;
double t85 = t82 + t83 + t84;
double t86 = 2*t20*t30*t85;
double t87 = 2*r*t20*t25*t33;
double t89 = -2*t25*t7*t88;
double t90 = t86 + t87 + t89;
double t91 = -2*t13*t54*t79*t80*t90;
double t92 = t78 + t91;
double t93 = sqrt(t92);
double t94 = t74 + t93;
double t2 = 1/sin(th);
double t3 = pow(t2,2);
gg = sqrt(-(((-t7 - t18*t8)*t88)/pow(t15,2)) + (pow(r,13)*spin*t16*t20*t55*t60*t94)/pow(t25,2) - (pow(r,20)*pow(t15,4)*t20*t3*t46*pow(t94,2))/(4.*pow(t54,2)))/(1 - (ktkp*pow(r,10)*t16*t3*t55*t94)/2.);
}
