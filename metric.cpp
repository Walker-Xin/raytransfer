void metric(double r, double chi, double g[4][4])
{
double t1 = pow(r,2);
double t2 = pow(chi,2);
double t3 = pow(spin,2);
double t4 = t2*t3;
double t5 = t1 + t4;
double t6 = 1/t5;
double t11 = 1/r;
double t12 = (defpar*spin*t11)/2.;
double t7 = -2*r;
double t13 = -1 + t2;
double t14 = 2*r*spin*t13*t6;
g[0][0] = -(t6*(t1 + t4 + t7));
g[0][1] = t12;
g[0][2] = 0;
g[0][3] = t14;
g[1][0] = t12;
g[1][1] = t5/(t1 + t3 + t7);
g[1][2] = 0;
g[1][3] = 0;
g[2][0] = 0;
g[2][1] = 0;
g[2][2] = t5/(1 - t2);
g[2][3] = 0;
g[3][0] = t14;
g[3][1] = 0;
g[3][2] = 0;
g[3][3] = -(t13*(pow(r,4) + pow(spin,4)*t2 + r*(2 + r - 2*t2 + r*t2)*t3)*t6);
}

void uppermetric(double r, double chi, double gu[4][4])
{
double t1 = pow(r,2);
double t4 = pow(spin,2);
double t7 = pow(chi,2);
double t11 = pow(r,4);
double t22 = pow(defpar,2);
double t12 = -2*t7;
double t13 = r*t7;
double t14 = 2 + r + t12 + t13;
double t16 = pow(spin,4);
double t15 = r*t14*t4;
double t17 = t16*t7;
double t18 = t11 + t15 + t17;
double t19 = pow(r,6);
double t20 = 4*t19;
double t21 = 8*t7;
double t23 = t21 + t22;
double t24 = t11*t23*t4;
double t25 = pow(chi,4);
double t26 = 4*r*t25;
double t27 = t14*t22;
double t28 = t26 + t27;
double t29 = r*t16*t28;
double t30 = pow(spin,6);
double t31 = t22*t30*t7;
double t32 = t20 + t24 + t29 + t31;
double t33 = 1/t32;
double t2 = -2 + r;
double t3 = r*t2;
double t5 = t3 + t4;
double t6 = 1/t5;
double t8 = t4*t7;
double t9 = t1 + t8;
double t35 = 2*defpar*r*spin*t18*t33;
double t36 = pow(r,3);
double t37 = -8*spin*t33*t36*t6*t9;
double t39 = 4*defpar*t1*t33*t4;
gu[0][0] = -4*t1*t18*t33*t6*t9;
gu[0][1] = t35;
gu[0][2] = 0;
gu[0][3] = t37;
gu[1][0] = t35;
gu[1][1] = 4*t1*t33*t5*t9;
gu[1][2] = 0;
gu[1][3] = t39;
gu[2][0] = 0;
gu[2][1] = 0;
gu[2][2] = (1 - t7)/t9;
gu[2][3] = 0;
gu[3][0] = t37;
gu[3][1] = t39;
gu[3][2] = 0;
gu[3][3] = -((t33*t6*(t16*t22 + 4*t2*t36 + r*t4*(t2*t22 + 4*r*t7))*t9)/(-1 + t7));
}

void metric_rderivatives(double r, double chi, double dg[4][4])
{
double t1 = pow(r,2);
double t3 = pow(chi,2);
double t4 = pow(spin,2);
double t5 = t3*t4;
double t2 = -t1;
double t6 = t2 + t5;
double t7 = t1 + t5;
double t8 = pow(t7,-2);
double t11 = pow(r,-2);
double t12 = -0.5*(defpar*spin*t11);
double t13 = -1 + t3;
double t14 = 2*spin*t13*t6*t8;
dg[0][0] = 2*t6*t8;
dg[0][1] = t12;
dg[0][2] = 0;
dg[0][3] = t14;
dg[1][0] = t12;
dg[1][1] = (-2*t1 + 2*(r + (1 - r)*t3)*t4)/pow((-2 + r)*r + t4,2);
dg[1][2] = 0;
dg[1][3] = 0;
dg[2][0] = 0;
dg[2][1] = 0;
dg[2][2] = (-2*r)/t13;
dg[2][3] = 0;
dg[3][0] = t14;
dg[3][1] = 0;
dg[3][2] = 0;
dg[3][3] = -2*t13*(pow(r,5) + pow(spin,4)*(pow(chi,4)*(-1 + r) + t3) + t1*(-1 + t3 + 2*r*t3)*t4)*t8;
}

void metric_r2derivatives(double r, double chi, double dg2[4][4])
{
double t1 = pow(r,2);
double t2 = pow(chi,2);
double t3 = pow(spin,2);
double t6 = t2*t3;
double t7 = t1 + t6;
double t8 = pow(t7,-3);
double t11 = pow(r,-3);
double t12 = defpar*spin*t11;
double t13 = -1 + t2;
double t14 = -t1;
double t15 = 3*t2*t3;
double t16 = t14 + t15;
double t17 = -4*r*spin*t13*t16*t8;
double t22 = pow(r,3);
double t31 = pow(spin,4);
dg2[0][0] = 4*r*(t1 - 3*t2*t3)*t8;
dg2[0][1] = t12;
dg2[0][2] = 0;
dg2[0][3] = t17;
dg2[1][0] = t12;
dg2[1][1] = (4*t22 + 2*(-3*t1 + (4 - 6*r + 3*t1)*t2)*t3 - 2*t13*t31)/pow((-2 + r)*r + t3,3);
dg2[1][2] = 0;
dg2[1][3] = 0;
dg2[2][0] = 0;
dg2[2][1] = 0;
dg2[2][2] = -2/t13;
dg2[2][3] = 0;
dg2[3][0] = t17;
dg2[3][1] = 0;
dg2[3][2] = 0;
dg2[3][3] = -2*t13*(pow(r,6) + pow(chi,6)*pow(spin,6) + (-2*t13 + 3*r*t2)*t22*t3 + 3*r*t2*(2*t13 + r*t2)*t31)*t8;
}