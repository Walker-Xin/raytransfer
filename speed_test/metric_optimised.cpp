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

	long double var1 = pow(sin(th),2);
    long double var2 = pow(spin,2);
    long double var3 = pow(r,2);
    long double var4 = pow(r,3);

    // gtt
	gtt = r*(var4 + r*var2*pow(cos(th),2))*pow((a13 + var4)*(var3 + var2) - r*(a22 + var3)*var2*var1,-2)*(-(pow(r,4)*((-2 + r)*r + var2)) + var2*pow(a22 + var3,2)*var1);
	// grr
    grr = r*(var4 + r*var2*pow(cos(th),2))*pow(a52 + var3,-1)*pow((-2 + r)*r + var2,-1);
	// gthth
    gthth = var3 + var2*pow(cos(th),2);
	// gpp
    gpp = pow(r,-1)*(var4 + r*var2*pow(cos(th),2))*pow((a13 + var4)*(var3 + var2) - r*(a22 + var3)*var2*var1,-2)*var1*(pow(a13 + var4,2)*pow(var3 + var2,2) - pow(r,6)*var2*((-2 + r)*r + var2)*var1);
	// gtp
    gtp = -2*spin*(2*r - var3 - var2 + pow(r,-5)*(a22 + var3)*(a13 + var4)*(var3 + var2))*(var3 + var2*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(var3 + var2) - (1 + a22*pow(r,-2))*var2*var1,-2)*var1;

	g[0][0] = gtt;
	g[0][3] = gtp;
	g[1][1] = grr;
	g[2][2] = gthth;
	g[3][0] = g[0][3];
	g[3][3] = gpp;
}

void metric_rderivatives(long double r, long double th, long double dg[][4])
{
	long double dgttdr, dgtpdr, dgppdr;

	long double var1 = pow(sin(th),2);
    long double var2 = pow(cos(th),2);
    long double var3 = pow(spin,2);
    long double var4 = pow(r,2);
    long double var5 = pow(r,3);

    // dgttdr
	dgttdr = r*(var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*(-((-2 + 2*r)*pow(r,4)) - 4*var5*((-2 + r)*r + var3) + 4*r*(a22 + var4)*var3*var1) + r*(3*var4 + var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*(-(pow(r,4)*((-2 + r)*r + var3)) + var3*pow(a22 + var4,2)*var1) + (var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*(-(pow(r,4)*((-2 + r)*r + var3)) + var3*pow(a22 + var4,2)*var1) - 2*r*(var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-3)*(2*r*(a13 + var5) + 3*var4*(var4 + var3) - 2*var4*var3*var1 - (a22 + var4)*var3*var1)*(-(pow(r,4)*((-2 + r)*r + var3)) + var3*pow(a22 + var4,2)*var1);
	// dgppdr
	dgppdr = pow(r,-1)*(var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*var1*(4*r*(var4 + var3)*pow(a13 + var5,2) + 6*var4*(a13 + var5)*pow(var4 + var3,2) - (-2 + 2*r)*pow(r,6)*var3*var1 - 6*pow(r,5)*var3*((-2 + r)*r + var3)*var1) + pow(r,-1)*(3*var4 + var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*var1*(pow(a13 + var5,2)*pow(var4 + var3,2) - pow(r,6)*var3*((-2 + r)*r + var3)*var1) - pow(r,-2)*(var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-2)*var1*(pow(a13 + var5,2)*pow(var4 + var3,2) - pow(r,6)*var3*((-2 + r)*r + var3)*var1) - 2*pow(r,-1)*(var5 + r*var3*var2)*pow((a13 + var5)*(var4 + var3) - r*(a22 + var4)*var3*var1,-3)*var1*(2*r*(a13 + var5) + 3*var4*(var4 + var3) - 2*var4*var3*var1 - (a22 + var4)*var3*var1)*(pow(a13 + var5,2)*pow(var4 + var3,2) - pow(r,6)*var3*((-2 + r)*r + var3)*var1);
	// dgtpdr
	dgtpdr = -4*r*spin*(2*r - var4 - var3 + pow(r,-5)*(a22 + var4)*(a13 + var5)*(var4 + var3))*pow((1 + a13*pow(r,-3))*(var4 + var3) - (1 + a22*pow(r,-2))*var3*var1,-2)*var1 - 2*spin*(2 - 2*r + 2*pow(r,-4)*(a22 + var4)*(a13 + var5) + 3*pow(r,-3)*(a22 + var4)*(var4 + var3) + 2*pow(r,-4)*(a13 + var5)*(var4 + var3) - 5*pow(r,-6)*(a22 + var4)*(a13 + var5)*(var4 + var3))*(var4 + var3*var2)*pow((1 + a13*pow(r,-3))*(var4 + var3) - (1 + a22*pow(r,-2))*var3*var1,-2)*var1 + 4*spin*(2*r - var4 - var3 + pow(r,-5)*(a22 + var4)*(a13 + var5)*(var4 + var3))*(var4 + var3*var2)*pow((1 + a13*pow(r,-3))*(var4 + var3) - (1 + a22*pow(r,-2))*var3*var1,-3)*var1*(2*r*(1 + a13*pow(r,-3)) - 3*a13*pow(r,-4)*(var4 + var3) + 2*a22*pow(r,-3)*var3*var1);

 	dg[0][0] = dgttdr;
	dg[0][3] = dgtpdr;
	dg[3][0] = dg[0][3];
	dg[3][3] = dgppdr;
}

void metric_r2derivatives(long double r, long double th, long double dg2[][4])
{
    long double dgttdr2, dgtpdr2, dgppdr2;

	long double v11 = pow(r,3);
	long double v12 = pow(spin,2);
	long double v13 = pow(cos(th),2);
	long double v14 = pow(r,2);
	long double v15 = pow(sin(th),2);
	long double v16 = pow(r,4);
	long double v17 = pow(r,-5);
	long double v18 = pow(r,-3);
	long double v19 = pow(r,-4);
	long double v110 = pow(r,-2);
	long double v111 = pow(r,-6);
	long double v112 = pow(r,-1);
	long double v113 = pow(r,5);
	long double v114 = pow(r,6);
	long double v21 = pow((a13 + v11)*(v14 + v12) - r*(a22 + v14)*v12*v15,-2);
	long double v22 = pow((a13 + v11)*(v14 + v12) - r*(a22 + v14)*v12*v15,-3);
	long double v23 = pow(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15,2);
	long double v24 = pow((a13 + v11)*(v14 + v12) - r*(a22 + v14)*v12*v15,-4);
	long double v25 = pow(a22 + v14,2);
	long double v26 = pow((1 + a13*v18)*(v14 + v12) - (1 + a22*v110)*v12*v15,-2);
	long double v27 = pow((1 + a13*v18)*(v14 + v12) - (1 + a22*v110)*v12*v15,-3);
	long double v28 = pow(a13 + v11,2);
	long double v29 = pow(v14 + v12,2);

    // dgttdr2
	dgttdr2 = r*(v11 + r*v12*v13)*v21*(-8*(-2 + 2*r)*v11 - 2*v16 - 12*v14*((-2 + r)*r + v12) + 8*v14*v12*v15 + 4*(a22 + v14)*v12*v15) + 2*r*(3*v14 + v12*v13)*v21*(-((-2 + 2*r)*v16) - 4*v11*((-2 + r)*r + v12) + 4*r*(a22 + v14)*v12*v15) + 2*(v11 + r*v12*v13)*v21*(-((-2 + 2*r)*v16) - 4*v11*((-2 + r)*r + v12) + 4*r*(a22 + v14)*v12*v15) - 4*r*(v11 + r*v12*v13)*v22*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(-((-2 + 2*r)*v16) - 4*v11*((-2 + r)*r + v12) + 4*r*(a22 + v14)*v12*v15) + 6*r*(v11 + r*v12*v13)*v23*v24*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15) + 6*v14*v21*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15) + 2*(3*v14 + v12*v13)*v21*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15) - 2*r*(v11 + r*v12*v13)*v22*(12*v11 + 2*(a13 + v11) + 6*r*(v14 + v12) - 6*r*v12*v15)*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15) - 4*r*(3*v14 + v12*v13)*v22*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15) - 4*(v11 + r*v12*v13)*v22*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(-(v16*((-2 + r)*r + v12)) + v12*v25*v15);
	//dgtpdr2
	dgtpdr2 = -12*spin*(2*r - v14 - v12 + v17*(a22 + v14)*(a13 + v11)*(v14 + v12))*(v14 + v12*v13)*pow(2*r*(1 + a13*v18) - 3*a13*v19*(v14 + v12) + 2*a22*v18*v12*v15,2)*pow((1 + a13*v18)*(v14 + v12) - (1 + a22*v110)*v12*v15,-4)*v15 - 8*r*spin*(2 - 2*r + 2*v19*(a22 + v14)*(a13 + v11) + 3*v18*(a22 + v14)*(v14 + v12) + 2*v19*(a13 + v11)*(v14 + v12) - 5*v111*(a22 + v14)*(a13 + v11)*(v14 + v12))*v26*v15 - 4*spin*(2*r - v14 - v12 + v17*(a22 + v14)*(a13 + v11)*(v14 + v12))*v26*v15 - 2*spin*(-2 + 12*v110*(a22 + v14) + 8*v18*(a13 + v11) - 18*v17*(a22 + v14)*(a13 + v11) + 12*v110*(v14 + v12) - 24*v19*(a22 + v14)*(v14 + v12) - 18*v17*(a13 + v11)*(v14 + v12) + 30*pow(r,-7)*(a22 + v14)*(a13 + v11)*(v14 + v12))*(v14 + v12*v13)*v26*v15 + 4*spin*(2*r - v14 - v12 + v17*(a22 + v14)*(a13 + v11)*(v14 + v12))*(v14 + v12*v13)*v27*v15*(-12*a13*v18 + 2*(1 + a13*v18) + 12*a13*v17*(v14 + v12) - 6*a22*v19*v12*v15) + 16*r*spin*(2*r - v14 - v12 + v17*(a22 + v14)*(a13 + v11)*(v14 + v12))*v27*v15*(2*r*(1 + a13*v18) - 3*a13*v19*(v14 + v12) + 2*a22*v18*v12*v15) + 8*spin*(2 - 2*r + 2*v19*(a22 + v14)*(a13 + v11) + 3*v18*(a22 + v14)*(v14 + v12) + 2*v19*(a13 + v11)*(v14 + v12) - 5*v111*(a22 + v14)*(a13 + v11)*(v14 + v12))*(v14 + v12*v13)*v27*v15*(2*r*(1 + a13*v18) - 3*a13*v19*(v14 + v12) + 2*a22*v18*v12*v15);
	// dgppdr2
	dgppdr2 = v112*(v11 + r*v12*v13)*v21*v15*(48*v11*(a13 + v11)*(v14 + v12) + 8*v14*v28 + 4*(v14 + v12)*v28 + 12*r*(a13 + v11)*v29 + 18*v16*v29 - 12*(-2 + 2*r)*v113*v12*v15 - 2*v114*v12*v15 - 30*v16*v12*((-2 + r)*r + v12)*v15) + 2*v112*(3*v14 + v12*v13)*v21*v15*(4*r*(v14 + v12)*v28 + 6*v14*(a13 + v11)*v29 - (-2 + 2*r)*v114*v12*v15 - 6*v113*v12*((-2 + r)*r + v12)*v15) - 2*v110*(v11 + r*v12*v13)*v21*v15*(4*r*(v14 + v12)*v28 + 6*v14*(a13 + v11)*v29 - (-2 + 2*r)*v114*v12*v15 - 6*v113*v12*((-2 + r)*r + v12)*v15) - 4*v112*(v11 + r*v12*v13)*v22*v15*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(4*r*(v14 + v12)*v28 + 6*v14*(a13 + v11)*v29 - (-2 + 2*r)*v114*v12*v15 - 6*v113*v12*((-2 + r)*r + v12)*v15) + 6*v112*(v11 + r*v12*v13)*v23*v24*v15*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) + 6*v21*v15*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) - 2*v110*(3*v14 + v12*v13)*v21*v15*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) + 2*v18*(v11 + r*v12*v13)*v21*v15*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) - 2*v112*(v11 + r*v12*v13)*v22*v15*(12*v11 + 2*(a13 + v11) + 6*r*(v14 + v12) - 6*r*v12*v15)*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) - 4*v112*(3*v14 + v12*v13)*v22*v15*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15) + 4*v110*(v11 + r*v12*v13)*v22*v15*(2*r*(a13 + v11) + 3*v14*(v14 + v12) - 2*v14*v12*v15 - (a22 + v14)*v12*v15)*(v28*v29 - v114*v12*((-2 + r)*r + v12)*v15);

	dg2[0][0] = dgttdr2;
	dg2[0][3] = dgtpdr2;
	dg2[3][0] = dg2[0][3];
	dg2[3][3] = dgppdr2;
}

void uppermetric(long double r, long double th, long double rth[])
{
	long double gurr, guthth;

	gurr = pow(r,-1)*(a52 + pow(r,2))*((-2 + r)*r + pow(spin,2))*pow(pow(r,3) + r*pow(spin,2)*pow(cos(th),2),-1);
	guthth = pow(pow(r,2) + pow(spin,2)*pow(cos(th),2),-1);

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