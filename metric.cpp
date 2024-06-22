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
double t15 = 2*r*spin*t13*t6;
g[0][0] = -(t6*(t1 + t4 + t7));
g[0][1] = t12;
g[0][2] = 0;
g[0][3] = t15;
g[1][0] = t12;
g[1][1] = t5/(t1 + t3 + t7);
g[1][2] = 0;
g[1][3] = 0;
g[2][0] = 0;
g[2][1] = 0;
g[2][2] = t5/(1 - t2);
g[2][3] = 0;
g[3][0] = t15;
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
double t23 = pow(defpar,2);
double t12 = -2*t7;
double t13 = r*t7;
double t15 = 2 + r + t12 + t13;
double t17 = pow(spin,4);
double t16 = r*t15*t4;
double t18 = t17*t7;
double t19 = t11 + t16 + t18;
double t20 = pow(r,6);
double t21 = 4*t20;
double t22 = 8*t7;
double t24 = t22 + t23;
double t25 = t11*t24*t4;
double t26 = pow(chi,4);
double t27 = 4*r*t26;
double t28 = t15*t23;
double t29 = t27 + t28;
double t30 = r*t17*t29;
double t31 = pow(spin,6);
double t32 = t23*t31*t7;
double t33 = t21 + t25 + t30 + t32;
double t34 = 1/t33;
double t2 = -2 + r;
double t3 = r*t2;
double t5 = t3 + t4;
double t6 = 1/t5;
double t8 = t4*t7;
double t9 = t1 + t8;
double t36 = 2*defpar*r*spin*t19*t34;
double t37 = pow(r,3);
double t38 = -8*spin*t34*t37*t6*t9;
double t40 = 4*defpar*t1*t34*t4;
gu[0][0] = -4*t1*t19*t34*t6*t9;
gu[0][1] = t36;
gu[0][2] = 0;
gu[0][3] = t38;
gu[1][0] = t36;
gu[1][1] = 4*t1*t34*t5*t9;
gu[1][2] = 0;
gu[1][3] = t40;
gu[2][0] = 0;
gu[2][1] = 0;
gu[2][2] = (1 - t7)/t9;
gu[2][3] = 0;
gu[3][0] = t38;
gu[3][1] = t40;
gu[3][2] = 0;
gu[3][3] = -((t34*t6*(t17*t23 + 4*t2*t37 + r*t4*(t2*t23 + 4*r*t7))*t9)/(-1 + t7));
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
double t15 = 2*spin*t13*t6*t8;
dg[0][0] = 2*t6*t8;
dg[0][1] = t12;
dg[0][2] = 0;
dg[0][3] = t15;
dg[1][0] = t12;
dg[1][1] = (-2*t1 + 2*(r + (1 - r)*t3)*t4)/pow((-2 + r)*r + t4,2);
dg[1][2] = 0;
dg[1][3] = 0;
dg[2][0] = 0;
dg[2][1] = 0;
dg[2][2] = (-2*r)/t13;
dg[2][3] = 0;
dg[3][0] = t15;
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
double t15 = -t1;
double t16 = 3*t2*t3;
double t17 = t15 + t16;
double t18 = -4*r*spin*t13*t17*t8;
double t23 = pow(r,3);
double t32 = pow(spin,4);
dg2[0][0] = 4*r*(t1 - 3*t2*t3)*t8;
dg2[0][1] = t12;
dg2[0][2] = 0;
dg2[0][3] = t18;
dg2[1][0] = t12;
dg2[1][1] = (4*t23 + 2*(-3*t1 + (4 - 6*r + 3*t1)*t2)*t3 - 2*t13*t32)/pow((-2 + r)*r + t3,3);
dg2[1][2] = 0;
dg2[1][3] = 0;
dg2[2][0] = 0;
dg2[2][1] = 0;
dg2[2][2] = -2/t13;
dg2[2][3] = 0;
dg2[3][0] = t18;
dg2[3][1] = 0;
dg2[3][2] = 0;
dg2[3][3] = -2*t13*(pow(r,6) + pow(chi,6)*pow(spin,6) + (-2*t13 + 3*r*t2)*t23*t3 + 3*r*t2*(2*t13 + r*t2)*t32)*t8;
}