void metric(double r, double th, double g[4][4])
{
double t3 = pow(spin,2);
double t4 = 2*th;
double t5 = cos(t4);
double t6 = t3*t5;
double t8 = pow(r,2);
double t9 = 2*t8;
double t12 = t3 + t6 + t9;
double t13 = 1/t12;
double t15 = 1/r;
double t16 = (defpar*spin*t15)/2.;
double t23 = cos(th);
double t24 = pow(t23,2);
double t25 = t24*t3;
double t26 = t25 + t8;
double t17 = sin(th);
double t18 = pow(t17,2);
double t19 = -4*r*spin*t13*t18;
double t1 = -2 + r;
g[0][0] = -(t13*(2*r*t1 + t3 + t6));
g[0][1] = t16;
g[0][2] = 0;
g[0][3] = t19;
g[1][0] = t16;
g[1][1] = t26/(-2*r + t3 + t8);
g[1][2] = 0;
g[1][3] = 0;
g[2][0] = 0;
g[2][1] = 0;
g[2][2] = t26;
g[2][3] = 0;
g[3][0] = t19;
g[3][1] = 0;
g[3][2] = 0;
g[3][3] = (t18*t26*(-(t18*t3*(r*t1 + t3)) + pow(t3 + t8,2)))/pow(t3 - t18*t3 + t8,2);
}

void uppermetric(double r, double th, double gu[4][4])
{
double t1 = pow(r,2);
double t4 = pow(spin,2);
double t2 = -2 + r;
double t3 = r*t2;
double t5 = t3 + t4;
double t8 = 2*th;
double t9 = cos(t8);
double t14 = pow(r,4);
double t25 = pow(defpar,2);
double t19 = pow(spin,4);
double t15 = 2*t14;
double t16 = 3*r;
double t17 = 2 + t16;
double t18 = r*t17*t4;
double t20 = t4*t5*t9;
double t21 = t15 + t18 + t19 + t20;
double t22 = pow(r,6);
double t23 = 8*t22;
double t24 = 8*t14*t4;
double t26 = 2*t14*t25*t4;
double t27 = 2*r*t19*t25;
double t28 = 3*t1*t19;
double t29 = 3*t1*t19*t25;
double t30 = pow(spin,6);
double t31 = t25*t30;
double t32 = 8*t14;
double t33 = t2*t25;
double t34 = 4*r;
double t35 = t33 + t34;
double t36 = r*t35*t4;
double t37 = t19*t25;
double t38 = t32 + t36 + t37;
double t39 = t38*t4*t9;
double t40 = 4*th;
double t41 = cos(t40);
double t42 = t1*t19*t41;
double t43 = t23 + t24 + t26 + t27 + t28 + t29 + t31 + t39 + t42;
double t44 = 1/t43;
double t6 = 1/t5;
double t7 = 2*t1;
double t12 = t4*t9;
double t13 = t12 + t4 + t7;
double t46 = 2*defpar*r*spin*t21*t44;
double t47 = pow(r,3);
double t48 = -8*spin*t13*t44*t47*t6;
double t50 = 8*defpar*t1*t4*t44;
gu[0][0] = -2*t1*t13*t21*t44*t6;
gu[0][1] = t46;
gu[0][2] = 0;
gu[0][3] = t48;
gu[1][0] = t46;
gu[1][1] = 4*t1*t13*t44*t5;
gu[1][2] = 0;
gu[1][3] = t50;
gu[2][0] = 0;
gu[2][1] = 0;
gu[2][2] = 1/(t1 + t4*pow(cos(th),2));
gu[2][3] = 0;
gu[3][0] = t48;
gu[3][1] = t50;
gu[3][2] = 0;
gu[3][3] = t13*t44*t6*(t37 + r*(2*r + t33)*t4 + 4*t2*t47 + 2*t1*t4*t9)*pow(1/sin(th),2);
}

void metric_rderivatives(double r, double th, double dg[4][4])
{
double t3 = pow(spin,2);
double t1 = pow(r,2);
double t4 = 2*th;
double t5 = cos(t4);
double t6 = t3*t5;
double t2 = -2*t1;
double t7 = t2 + t3 + t6;
double t8 = 2*t1;
double t9 = t3 + t6 + t8;
double t12 = pow(t9,-2);
double t14 = pow(r,-2);
double t15 = -0.5*(defpar*spin*t14);
double t23 = -r;
double t16 = sin(th);
double t17 = pow(t16,2);
double t18 = -4*spin*t12*t17*t7;
double t38 = pow(spin,4);
double t32 = 2*r;
double t26 = 1 + t23;
dg[0][0] = 4*t12*t7;
dg[0][1] = t15;
dg[0][2] = 0;
dg[0][3] = t18;
dg[1][0] = t15;
dg[1][1] = (2*(r*(t23 + t3) + t26*t3*pow(cos(th),2)))/pow((-2 + r)*r + t3,2);
dg[1][2] = 0;
dg[1][3] = 0;
dg[2][0] = 0;
dg[2][1] = 0;
dg[2][2] = t32;
dg[2][3] = 0;
dg[3][0] = t18;
dg[3][1] = 0;
dg[3][2] = 0;
dg[3][3] = t12*t17*(8*pow(r,5) + 8*pow(r,3)*t3 - 4*t1*t3 + t38 + 3*r*t38 + 4*r*t3*(t3 + r*(1 + t32))*t5 - t26*t38*cos(4*th));
}

void metric_r2derivatives(double r, double th, double dg2[4][4])
{
double t3 = pow(spin,2);
double t1 = pow(r,2);
double t2 = 2*t1;
double t5 = 2*th;
double t6 = cos(t5);
double t9 = t3*t6;
double t12 = t2 + t3 + t9;
double t13 = pow(t12,-3);
double t15 = pow(r,-3);
double t16 = defpar*spin*t15;
double t17 = -2*t1;
double t18 = 3*t3;
double t19 = 3*t3*t6;
double t20 = t17 + t18 + t19;
double t21 = sin(th);
double t22 = pow(t21,2);
double t23 = 16*r*spin*t13*t20*t22;
double t28 = pow(r,3);
double t31 = pow(spin,4);
double t47 = pow(spin,6);
dg2[0][0] = 16*r*t13*(t2 - 3*t3 - 3*t3*t6);
dg2[0][1] = t16;
dg2[0][2] = 0;
dg2[0][3] = t23;
dg2[1][0] = t16;
dg2[1][1] = (2*(2*t28 - 3*t1*t3 + t31 - t3*(-4 + 6*r - 3*t1 + t3)*pow(cos(th),2)))/pow((-2 + r)*r + t3,3);
dg2[1][2] = 0;
dg2[1][3] = 0;
dg2[2][0] = 0;
dg2[2][1] = 0;
dg2[2][2] = 2;
dg2[2][3] = 0;
dg2[3][0] = t23;
dg2[3][1] = 0;
dg2[3][2] = 0;
dg2[3][3] = (t13*t22*(32*pow(r,6) + 48*pow(r,4)*t3 + 32*t28*t3 - 24*r*t31 + 36*t1*t31 + 10*t47 + t3*(16*(-2 + 3*r)*t28 + 48*t1*t3 + 15*t31)*t6 + 6*(2*r*(2 + r) + t3)*t31*cos(4*th) + t47*cos(6*th)))/2.;
}