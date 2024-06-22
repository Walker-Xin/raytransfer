void redshift(double r, double chi, double ktkp, double &gg)
{
double t1 = pow(chi,2);
double t4 = pow(r,2);
double t5 = pow(spin,2);
double t2 = -1 + t1;
double t6 = t1*t5;
double t7 = t4 + t6;
double t23 = -t4;
double t24 = t23 + t6;
double t29 = pow(t7,-4);
double t9 = pow(r,5);
double t11 = 2*r*t1;
double t12 = -1 + t1 + t11;
double t13 = t12*t4*t5;
double t15 = pow(chi,4);
double t16 = -1 + r;
double t17 = t15*t16;
double t18 = t1 + t17;
double t19 = pow(spin,4);
double t20 = t18*t19;
double t21 = t13 + t20 + t9;
double t22 = 1/t21;
double t25 = pow(t7,-2);
double t26 = -2*spin*t2*t24*t25;
double t27 = pow(t2,2);
double t28 = pow(t24,2);
double t30 = 4*t27*t28*t29*t5;
double t31 = 4*t2*t21*t24*t29;
double t32 = t30 + t31;
double t33 = sqrt(t32);
double t34 = t26 + t33;
double t3 = 1/t2;
gg = sqrt((-2*r + t4 + t6)/t7 + 2*r*spin*t22*t34*t7 + (t3*pow(t34,2)*(pow(r,4) + t1*t19 + r*(2 + r - 2*t1 + r*t1)*t5)*pow(t7,3))/(4.*pow(t21,2)))/(1 + (ktkp*t22*t3*t34*pow(t7,2))/2.);
}

double specific_energy(double r)
{
double t3 = pow(r,2);
double t1 = pow(r,-2);
double t8 = pow(spin,2);
double t7 = pow(r,5);
double t9 = -(t3*t8);
double t11 = t7 + t9;
double t2 = -2*r;
double t4 = t2 + t3;
double t5 = t1*t4;
double t6 = pow(r,3);
double t12 = 1/t11;
double t13 = -2*spin*t1;
double t15 = pow(r,-4);
double t16 = 4*t15*t8;
double t17 = pow(r,-6);
double t18 = 4*t11*t17;
double t19 = t16 + t18;
double t20 = sqrt(t19);
double t21 = t13 + t20;
double se = (t5 + spin*t12*t21*t6)/sqrt(t5 + 2*spin*t12*t21*t6 - (pow(r,6)*pow(t21,2)*(pow(r,4) + r*(2 + r)*t8))/(4.*pow(t11,2)));
return se;
}

double specific_momentum(double r)
{
double t3 = pow(r,2);
double t5 = pow(spin,2);
double t4 = pow(r,5);
double t6 = -(t3*t5);
double t7 = t4 + t6;
double t15 = pow(r,-2);
double t8 = 1/t7;
double t16 = -2*spin*t15;
double t17 = pow(r,-4);
double t18 = 4*t17*t5;
double t19 = pow(r,-6);
double t20 = 4*t19*t7;
double t21 = t18 + t20;
double t22 = sqrt(t21);
double t23 = t16 + t22;
double t9 = pow(r,4);
double t11 = 2 + r;
double t12 = r*t11*t5;
double t13 = t12 + t9;
double Lz = ((2*spin)/r - (t13*t23*t3*t8)/2.)/sqrt(t15*(-2*r + t3) - (pow(r,6)*t13*pow(t23,2))/(4.*pow(t7,2)) + 2*pow(r,3)*spin*t23*t8);
return Lz;
}