#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <time.h>
using namespace std;

#define imax 90 //Number of radial values

/* RK45 constants */
#define a1 1.0/4.0
#define b1_rk 3.0/32.0
#define b2_rk 9.0/32.0
#define c1_rk 1932.0/2197.0
#define c2_rk -7200.0/2197.0
#define c3 7296.0/2197.0
#define d1 439.0/216.0
#define d2 -8.0
#define d3 3680.0/513.0
#define d4 -845.0/4104.0
#define e1 -8.0/27.0
#define e2 2.0
#define e3 -3544.0/2565.0
#define e4 1859.0/4104.0
#define e5 -11.0/40.0
#define f1 25.0/216.0
#define f2 0.0
#define f3 1408.0/2565.0
#define f4 2197.0/4104.0
#define f5 -1.0/5.0
#define g1 16.0/135.0
#define g2 0.0
#define g3 6656.0/12825.0
#define g4 28561.0/56430.0
#define g5 -9.0/50.0
#define g6 2.0/55.0

const double Pi = 3.14159265358979;

double defpar;
double spin;
double spin2 = spin*spin;
double isco;
double dobs, iobs_deg, inc;
double Mdl, eta;

void christoffel(double r, double th, double christ[4][4][4]);
void diffeqs(double vars[], double diffs[]);
void redshift(double r, double th, double ktkp, double &gg);
double find_isco(double spin, double defpar);
void intersection(double x_1, double y_1, double z_1, double x_2, double y_2, double z_2, double x_d[]);
void metric(double r, double th, double g[4][4]);
void uppermetric(double r, double th, double gu[4][4]);
void metric_rderivatives(double r, double th, double dg[4][4]);
void metric_r2derivatives(double r, double th, double dg2[4][4]);
int raytrace(double errmin, double errmax, double xscr, double yscr, double traced[4]);
int raytrace_RKN(double errtol, double xobs, double yobs, double traced[4]);

#include "metric.cpp"
#include "redshift.cpp"
#include "christoffel.cpp"
#include "diffeqs.cpp"
#include "find_isco.cpp"
#include "intersection.cpp"
#include "raytracing.cpp"

#endif
