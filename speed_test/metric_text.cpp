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

	// gtt
	gtt = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2));
	// grr
    grr = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow(a52 + pow(r,2),-1)*pow((-2 + r)*r + pow(spin,2),-1);
	// gthth
    gthth = pow(r,2) + pow(spin,2)*pow(cos(th),2);
	// gpp
    gpp = pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2));
	// gtp
    gtp = -2*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2);

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

	// dgttdr
	dgttdr = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-((-2 + 2*r)*pow(r,4)) - 4*pow(r,3)*((-2 + r)*r + pow(spin,2)) + 4*r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2)) + r*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) + (pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) - 2*r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2));
	// dgppdr
	dgppdr = pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(4*r*(pow(r,2) + pow(spin,2))*pow(a13 + pow(r,3),2) + 6*pow(r,2)*(a13 + pow(r,3))*pow(pow(r,2) + pow(spin,2),2) - (-2 + 2*r)*pow(r,6)*pow(spin,2)*pow(sin(th),2) - 6*pow(r,5)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + pow(r,-1)*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - pow(r,-2)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 2*pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2));
	// dgtpdr
	dgtpdr = -4*r*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2) - 2*spin*(2 - 2*r + 2*pow(r,-4)*(a22 + pow(r,2))*(a13 + pow(r,3)) + 3*pow(r,-3)*(a22 + pow(r,2))*(pow(r,2) + pow(spin,2)) + 2*pow(r,-4)*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - 5*pow(r,-6)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2) + 4*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(1 + a13*pow(r,-3)) - 3*a13*pow(r,-4)*(pow(r,2) + pow(spin,2)) + 2*a22*pow(r,-3)*pow(spin,2)*pow(sin(th),2));

 	dg[0][0] = dgttdr;
	dg[0][3] = dgtpdr;
	dg[3][0] = dg[0][3];
	dg[3][3] = dgppdr;
}

void metric_r2derivatives(long double r, long double th, long double dg2[][4])
{
	long double dgttdr2, dgtpdr2, dgppdr2;

	// dgttdr2
	dgttdr2 = r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-8*(-2 + 2*r)*pow(r,3) - 2*pow(r,4) - 12*pow(r,2)*((-2 + r)*r + pow(spin,2)) + 8*pow(r,2)*pow(spin,2)*pow(sin(th),2) + 4*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2)) + 2*r*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-((-2 + 2*r)*pow(r,4)) - 4*pow(r,3)*((-2 + r)*r + pow(spin,2)) + 4*r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2)) + 2*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-((-2 + 2*r)*pow(r,4)) - 4*pow(r,3)*((-2 + r)*r + pow(spin,2)) + 4*r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2)) - 4*r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(-((-2 + 2*r)*pow(r,4)) - 4*pow(r,3)*((-2 + r)*r + pow(spin,2)) + 4*r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2)) + 6*r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),2)*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-4)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) + 6*pow(r,2)*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) + 2*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) - 2*r*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*(12*pow(r,3) + 2*(a13 + pow(r,3)) + 6*r*(pow(r,2) + pow(spin,2)) - 6*r*pow(spin,2)*pow(sin(th),2))*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) - 4*r*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2)) - 4*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(-(pow(r,4)*((-2 + r)*r + pow(spin,2))) + pow(spin,2)*pow(a22 + pow(r,2),2)*pow(sin(th),2));
	//dgtpdr2
	dgtpdr2 = -12*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow(2*r*(1 + a13*pow(r,-3)) - 3*a13*pow(r,-4)*(pow(r,2) + pow(spin,2)) + 2*a22*pow(r,-3)*pow(spin,2)*pow(sin(th),2),2)*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-4)*pow(sin(th),2) - 8*r*spin*(2 - 2*r + 2*pow(r,-4)*(a22 + pow(r,2))*(a13 + pow(r,3)) + 3*pow(r,-3)*(a22 + pow(r,2))*(pow(r,2) + pow(spin,2)) + 2*pow(r,-4)*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - 5*pow(r,-6)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2) - 4*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2) - 2*spin*(-2 + 12*pow(r,-2)*(a22 + pow(r,2)) + 8*pow(r,-3)*(a13 + pow(r,3)) - 18*pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3)) + 12*pow(r,-2)*(pow(r,2) + pow(spin,2)) - 24*pow(r,-4)*(a22 + pow(r,2))*(pow(r,2) + pow(spin,2)) - 18*pow(r,-5)*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) + 30*pow(r,-7)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2) + 4*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(-12*a13*pow(r,-3) + 2*(1 + a13*pow(r,-3)) + 12*a13*pow(r,-5)*(pow(r,2) + pow(spin,2)) - 6*a22*pow(r,-4)*pow(spin,2)*pow(sin(th),2)) + 16*r*spin*(2*r - pow(r,2) - pow(spin,2) + pow(r,-5)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(1 + a13*pow(r,-3)) - 3*a13*pow(r,-4)*(pow(r,2) + pow(spin,2)) + 2*a22*pow(r,-3)*pow(spin,2)*pow(sin(th),2)) + 8*spin*(2 - 2*r + 2*pow(r,-4)*(a22 + pow(r,2))*(a13 + pow(r,3)) + 3*pow(r,-3)*(a22 + pow(r,2))*(pow(r,2) + pow(spin,2)) + 2*pow(r,-4)*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - 5*pow(r,-6)*(a22 + pow(r,2))*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)))*(pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((1 + a13*pow(r,-3))*(pow(r,2) + pow(spin,2)) - (1 + a22*pow(r,-2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(1 + a13*pow(r,-3)) - 3*a13*pow(r,-4)*(pow(r,2) + pow(spin,2)) + 2*a22*pow(r,-3)*pow(spin,2)*pow(sin(th),2));
	// dgppdr2
	dgppdr2 = pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(48*pow(r,3)*(a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) + 8*pow(r,2)*pow(a13 + pow(r,3),2) + 4*(pow(r,2) + pow(spin,2))*pow(a13 + pow(r,3),2) + 12*r*(a13 + pow(r,3))*pow(pow(r,2) + pow(spin,2),2) + 18*pow(r,4)*pow(pow(r,2) + pow(spin,2),2) - 12*(-2 + 2*r)*pow(r,5)*pow(spin,2)*pow(sin(th),2) - 2*pow(r,6)*pow(spin,2)*pow(sin(th),2) - 30*pow(r,4)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + 2*pow(r,-1)*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(4*r*(pow(r,2) + pow(spin,2))*pow(a13 + pow(r,3),2) + 6*pow(r,2)*(a13 + pow(r,3))*pow(pow(r,2) + pow(spin,2),2) - (-2 + 2*r)*pow(r,6)*pow(spin,2)*pow(sin(th),2) - 6*pow(r,5)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 2*pow(r,-2)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(4*r*(pow(r,2) + pow(spin,2))*pow(a13 + pow(r,3),2) + 6*pow(r,2)*(a13 + pow(r,3))*pow(pow(r,2) + pow(spin,2),2) - (-2 + 2*r)*pow(r,6)*pow(spin,2)*pow(sin(th),2) - 6*pow(r,5)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 4*pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(4*r*(pow(r,2) + pow(spin,2))*pow(a13 + pow(r,3),2) + 6*pow(r,2)*(a13 + pow(r,3))*pow(pow(r,2) + pow(spin,2),2) - (-2 + 2*r)*pow(r,6)*pow(spin,2)*pow(sin(th),2) - 6*pow(r,5)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + 6*pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),2)*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-4)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + 6*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 2*pow(r,-2)*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + 2*pow(r,-3)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-2)*pow(sin(th),2)*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 2*pow(r,-1)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(12*pow(r,3) + 2*(a13 + pow(r,3)) + 6*r*(pow(r,2) + pow(spin,2)) - 6*r*pow(spin,2)*pow(sin(th),2))*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) - 4*pow(r,-1)*(3*pow(r,2) + pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2)) + 4*pow(r,-2)*(pow(r,3) + r*pow(spin,2)*pow(cos(th),2))*pow((a13 + pow(r,3))*(pow(r,2) + pow(spin,2)) - r*(a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2),-3)*pow(sin(th),2)*(2*r*(a13 + pow(r,3)) + 3*pow(r,2)*(pow(r,2) + pow(spin,2)) - 2*pow(r,2)*pow(spin,2)*pow(sin(th),2) - (a22 + pow(r,2))*pow(spin,2)*pow(sin(th),2))*(pow(a13 + pow(r,3),2)*pow(pow(r,2) + pow(spin,2),2) - pow(r,6)*pow(spin,2)*((-2 + r)*r + pow(spin,2))*pow(sin(th),2));

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
	// run metric for 100000 times
	for (int i = 0; i < 100000; i++)
	{
		metric_r2derivatives(1.5, 3.1415926/2., m);
	}
	end = time(NULL);
	dif = difftime(end, start);
	printf("time: %f\n", dif);
}