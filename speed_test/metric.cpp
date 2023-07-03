/* Metric and derivative of metric used throughout code */
/* If possible, use Maple or Mathematica to optimize code */

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>
using namespace std;

long double epsi3 = 0.0;
long double a13 = 0.01;
long double a22 = 0.0;
long double a52 = 0.0;
long double spin = 0.5;

void metric(long double r, long double th, long double g[][4])
{
	long double gtt, grr, gthth, gpp, gtp;

	long double t2 = pow(r,3);
	long double t6 = pow(spin,2);
	long double t15 = pow(r,2);
	long double t18 = a22 + t15;
	long double t19 = sin(th);
	long double t20 = pow(t19,2);
	long double t25 = -2 + r;
	long double t26 = r*t25;
	long double t33 = t26 + t6;
	long double t7 = cos(th);
	long double t8 = pow(t7,2);
	long double t9 = r*t6*t8;
	long double t10 = epsi3 + t2 + t9;
	long double t43 = 1/r;
	long double t13 = a13 + t2;
	long double t16 = t15 + t6;
	long double t17 = t13*t16;
	long double t21 = -(r*t18*t20*t6);
	long double t22 = t17 + t21;
	long double t23 = pow(t22,-2);
	long double t44 = epsi3*t43;
	long double t46 = t6*t8;
	long double t47 = t15 + t44 + t46;

	gtt = r*t10*t23*(-(pow(r,4)*t33) + pow(t18,2)*t20*t6);
	grr = (r*t10)/((a52 + t15)*t33);
	gthth = t47;
	gpp = t10*t20*t23*t43*(pow(t13,2)*pow(t16,2) - pow(r,6)*t20*t33*t6);
	gtp = -((spin*t20*t47*(2*r - t15 + (t13*t16*t18)/pow(r,5) - t6))/pow((1 + a13/pow(r,3))*t16 - (1 + a22/pow(r,2))*t20*t6,2));

	g[0][0] = gtt;
	g[0][3] = gtp;
	g[1][1] = grr;
	g[2][2] = gthth;
	g[3][0] = g[0][3];
	g[3][3] = gpp;
}

void metric_rderivatives(long double r, long double th, long double dg[][4])
{
	long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29;
	long double dgttdr, dgtpdr, dgppdr;

	t1 = r * r;
	t2 = pow(t1, 0.2e1);
	t3 = r * t1;
	t4 = t1 * t2;
	t5 = a22 + t1;
	t6 = sin(th);
	t6 = pow(t6, 0.2e1);
	t7 = spin * spin;
	t8 = t5 * t6;
	t9 = cos(th);
	t9 = t7 * pow(t9, 0.2e1);
	t10 = (t1 + t9) * r + epsi3;
	t11 = t1 + t7;
	t12 = a13 + t3;
	t13 = t11 * t12;
	t14 = t8 * t7 * r - t13;
	t15 = 2;
	t5 = (double) t15 * r * t2 - t7 * (-t6 * pow(t5, 0.2e1) + t2) - t4;
	t16 = 0.3e1 * t1;
	t17 = t16 + t9;
	t18 = (0.5e1 / 0.2e1 * t1 + 0.3e1 / 0.2e1 * t7) * r + a13;
	t16 = r * t18 - t7 * (a22 + t16) * t6 / 0.2e1;
	t14 = 0.1e1 / t14;
	t19 = pow(t14, 0.2e1);
	t4 = t6 * ((double) t15 * t3 * t2 * t7 - t4 * t7 * t11) + pow(t11, 0.2e1) * pow(t12, 0.2e1);
	t12 = 0.1e1 / r;
	t20 = t10 * t12;
	t21 = t7 * a22;
	t22 = 0.3e1 / 0.5e1;
	t23 = pow(t12, 0.2e1);
	t24 = t12 * t23;
	t25 = a13 * t24 + 0.1e1;
	t26 = a22 * t23 + 0.1e1;
	t27 = t11 * t25;
	t28 = (double) t15 * r;
	t29 = -t7 * t26 * t6 + t27;
	t29 = 0.1e1 / t29;

	dgttdr = 0.4e1 * t10 * r * t19 * (r * (t3 * (-0.3e1 / 0.2e1 * r + 0.5e1 / 0.2e1) - t7 * (t1 - t8)) + t16 * t5 * t14) + t5 * t19 * (r * t17 + t10);
	dgppdr = 0.4e1 * t6 * t10 * t19 * (t13 * t18 + (0.7e1 / 0.2e1 * (-0.4e1 / 0.7e1 * r + 0.1e1) * r - 0.3e1 / 0.2e1 * t7) * t7 * t2 * t6 + t4 * t16 * t14 * t12) + t6 * t4 * t19 * t12 * (t17 - t20);
	dgtpdr = t6 * spin * (0.5e1 * t20 * (a13 * (t1 * (t1 / 0.5e1 + (a22 + t7) * t22) + t21) - 0.2e1 / 0.5e1 * t3 * (t3 - t21)) * t19 + (t27 * t26 - t1 + t28 - t7) * pow(t29, 0.2e1) * (epsi3 * t23 - t28 + (double) t15 * (epsi3 * t12 + t1 + t9) * t29 * ((double) t15 * (t21 * t24 * t6 + r * t25) - 0.3e1 * t11 * a13 * pow(t23, 0.2e1))));

 	dg[0][0] = dgttdr;
	dg[0][3] = dgtpdr;
	dg[3][0] = dg[0][3];
	dg[3][3] = dgppdr;
}

void metric_r2derivatives(long double r, long double th, long double dg2[][4])
{
	long double dgttdr2, dgtpdr2, dgppdr2;

	long double t2 = pow(r,3);
	long double t6 = pow(spin,2);
	long double t19 = pow(r,2);
	long double t24 = sin(th);
	long double t25 = pow(t24,2);
	long double t33 = a22 + t19;
	long double t7 = cos(th);
	long double t8 = pow(t7,2);
	long double t9 = r*t6*t8;
	long double t10 = epsi3 + t2 + t9;
	long double t36 = a13 + t2;
	long double t37 = t19 + t6;
	long double t38 = t36*t37;
	long double t39 = -(r*t25*t33*t6);
	long double t40 = t38 + t39;
	long double t13 = pow(r,4);
	long double t16 = 2*r;
	long double t17 = -2 + t16;
	long double t20 = -2 + r;
	long double t21 = r*t20;
	long double t22 = t21 + t6;
	long double t41 = pow(t40,-2);
	long double t50 = -(t13*t17);
	long double t51 = -4*t2*t22;
	long double t52 = 4*r*t25*t33*t6;
	long double t53 = t50 + t51 + t52;
	long double t43 = 2*r*t36;
	long double t44 = 3*t19*t37;
	long double t46 = -2*t19*t25*t6;
	long double t47 = -(t25*t33*t6);
	long double t48 = t43 + t44 + t46 + t47;
	long double t49 = pow(t40,-3);
	long double t62 = -(t13*t22);
	long double t63 = pow(t33,2);
	long double t64 = t25*t6*t63;
	long double t67 = t62 + t64;
	long double t55 = 3*t19;
	long double t56 = t6*t8;
	long double t57 = t55 + t56;
	long double t85 = pow(t36,2);
	long double t89 = pow(t37,2);
	long double t84 = 1/r;
	long double t92 = pow(r,6);
	long double t94 = pow(r,5);
	long double t99 = 4*r*t37*t85;
	long double t100 = 6*t19*t36*t89;
	long double t101 = -(t17*t25*t6*t92);
	long double t102 = -6*t22*t25*t6*t94;
	long double t103 = t100 + t101 + t102 + t99;
	long double t60 = pow(t48,2);
	long double t61 = pow(t40,-4);
	long double t70 = 12*t2;
	long double t71 = 2*t36;
	long double t73 = 6*r*t37;
	long double t74 = -6*r*t25*t6;
	long double t75 = t70 + t71 + t73 + t74;
	long double t109 = t85*t89;
	long double t110 = -(t22*t25*t6*t92);
	long double t111 = t109 + t110;
	long double t107 = pow(r,-2);
	long double t119 = pow(r,-3);
	long double t125 = pow(r,-4);
	long double t132 = a13*t119;
	long double t133 = 1 + t132;
	long double t134 = t133*t37;
	long double t135 = a22*t107;
	long double t136 = 1 + t135;
	long double t137 = -(t136*t25*t6);
	long double t138 = t134 + t137;
	long double t139 = pow(t138,-2);
	long double t145 = pow(r,-5);
	long double t143 = -t19;
	long double t144 = -t6;
	long double t146 = t145*t33*t36*t37;
	long double t147 = t143 + t144 + t146 + t16;
	long double t158 = epsi3*t84;
	long double t159 = t158 + t19 + t56;
	long double t122 = -(epsi3*t107);
	long double t123 = t122 + t16;
	long double t161 = pow(t138,-3);
	long double t124 = -2*r;
	long double t126 = 2*t125*t33*t36;
	long double t127 = 3*t119*t33*t37;
	long double t128 = 2*t125*t36*t37;
	long double t129 = pow(r,-6);
	long double t130 = -5*t129*t33*t36*t37;
	long double t131 = 2 + t124 + t126 + t127 + t128 + t130;
	long double t168 = 2*r*t133;
	long double t169 = -3*a13*t125*t37;
	long double t170 = 2*a22*t119*t25*t6;
	long double t171 = t168 + t169 + t170;
	
	dgttdr2 = 2*t10*t41*t53 - 4*r*t10*t48*t49*t53 + 2*r*t41*t53*t57 + r*t10*t41*(-2*t13 - 8*t17*t2 - 12*t19*t22 + 8*t19*t25*t6 + 4*t25*t33*t6) + 6*t19*t41*t67 - 4*t10*t48*t49*t67 + 2*t41*t57*t67 - 4*r*t48*t49*t57*t67 + 6*r*t10*t60*t61*t67 - 2*r*t10*t49*t67*t75;
	dgppdr2 = -2*t10*t103*t107*t25*t41 + 6*t111*t25*t41 + 2*t10*t111*t119*t25*t41 + 4*t10*t107*t111*t25*t48*t49 - 2*t107*t111*t25*t41*t57 - 4*t10*t103*t25*t48*t49*t84 + 2*t103*t25*t41*t57*t84 - 4*t111*t25*t48*t49*t57*t84 + 6*t10*t111*t25*t60*t61*t84 - 2*t10*t111*t25*t49*t75*t84 + t10*t25*t41*t84*(48*t2*t36*t37 - 30*t13*t22*t25*t6 + 8*t19*t85 + 4*t37*t85 + 18*t13*t89 + 12*r*t36*t89 - 2*t25*t6*t92 - 12*t17*t25*t6*t94);
	dgtpdr2 = -4*spin*t123*t131*t139*t25 - 2*spin*(2 + 2*epsi3*t119)*t139*t147*t25 + 8*spin*t123*t147*t161*t171*t25 + 8*spin*t131*t159*t161*t171*t25 - (12*spin*t147*t159*pow(t171,2)*t25)/pow(t138,4) - 2*spin*t139*t159*t25*(-2 + 12*t107*t33 + 8*t119*t36 - 18*t145*t33*t36 + 12*t107*t37 - 24*t125*t33*t37 - 18*t145*t36*t37 + (30*t33*t36*t37)/pow(r,7)) + 4*spin*t147*t159*t161*t25*(-12*a13*t119 + 2*t133 + 12*a13*t145*t37 - 6*a22*t125*t25*t6);


	
	dg2[0][0] = dgttdr2;
	dg2[0][3] = dgtpdr2;
	dg2[3][0] = dg2[0][3];
	dg2[3][3] = dgppdr2;
}

void uppermetric(long double r, long double th, long double rth[])
{
	long double gurr, guthth;
	long double t1, t2, t3, t4;

	t1 = cos(th);
	t2 = 0.1e1 / r;
	t3 = r * r;
	t4 = spin * spin;
	t1 = t4 * pow(t1, 0.2e1) + epsi3 * t2 + t3;
	t1 = 0.1e1 / t1;

	gurr = t1 * (-0.2e1 * r + t3 + t4) * (a52 * pow(t2, 0.2e1) + 0.1e1);
	guthth = t1;

	rth[0] = gurr;
	rth[1] = guthth;
}

// testing computation time of metric
int main()
{
	time_t start, end;
	double dif;
	long double m[4][4];

	start = time(NULL);
	// run metric for 500000 times
	for (int i = 0; i < 500000; i++)
	{
		metric_r2derivatives(1.5, 3.1415926/2., m);
	}
	end = time(NULL);
	dif = difftime(end, start);
	printf("time: %f\n", dif);
}