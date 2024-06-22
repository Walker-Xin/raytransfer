void redshift(double r, double chi, double ktkp, double &gg)
{
double t1 = pow(chi,2);
double t4 = pow(r,2);
double t5 = pow(spin,2);
double t2 = -1 + t1;
double t6 = t1*t5;
double t7 = t4 + t6;
double t22 = -t4;
double t23 = t22 + t6;
double t28 = pow(t7,-4);
double t9 = pow(r,5);
double t11 = 2*r*t1;
double t12 = -1 + t1 + t11;
double t13 = t12*t4*t5;
double t14 = pow(chi,4);
double t15 = -1 + r;
double t16 = t14*t15;
double t17 = t1 + t16;
double t18 = pow(spin,4);
double t19 = t17*t18;
double t20 = t13 + t19 + t9;
double t21 = 1/t20;
double t24 = pow(t7,-2);
double t25 = -2*spin*t2*t23*t24;
double t26 = pow(t2,2);
double t27 = pow(t23,2);
double t29 = 4*t26*t27*t28*t5;
double t30 = 4*t2*t20*t23*t28;
double t31 = t29 + t30;
double t32 = sqrt(t31);
double t33 = t25 + t32;
double t3 = 1/t2;
gg = sqrt((-2*r + t4 + t6)/t7 + 2*r*spin*t21*t33*t7 + (t3*pow(t33,2)*(pow(r,4) + t1*t18 + r*(2 + r - 2*t1 + r*t1)*t5)*pow(t7,3))/(4.*pow(t20,2)))/(1 + (ktkp*t21*t3*t33*pow(t7,2))/2.);
}

double specific_energy(double r)
{
double t1 = pow(r,2);
double t2 = pow(chi,2);
double t3 = pow(spin,2);
double t4 = t2*t3;
double t5 = t1 + t4;
double t23 = -1 + t2;
double t24 = -t1;
double t25 = t24 + t4;
double t30 = pow(t5,-4);
double t11 = pow(r,5);
double t12 = 2*r*t2;
double t13 = -1 + t12 + t2;
double t14 = t1*t13*t3;
double t15 = pow(chi,4);
double t16 = -1 + r;
double t17 = t15*t16;
double t18 = t17 + t2;
double t19 = pow(spin,4);
double t20 = t18*t19;
double t21 = t11 + t14 + t20;
double t6 = 1/t5;
double t7 = -2*r;
double t8 = t1 + t4 + t7;
double t9 = t6*t8;
double t22 = 1/t21;
double t26 = pow(t5,-2);
double t27 = -2*spin*t23*t25*t26;
double t28 = pow(t23,2);
double t29 = pow(t25,2);
double t31 = 4*t28*t29*t3*t30;
double t32 = 4*t21*t23*t25*t30;
double t33 = t31 + t32;
double t34 = sqrt(t33);
double t35 = t27 + t34;
double se = (r*spin*t22*t35*t5 + t9)/sqrt(2*r*spin*t22*t35*t5 + ((pow(r,4) + t19*t2 + r*(2 + r - 2*t2 + r*t2)*t3)*pow(t35,2)*pow(t5,3))/(4.*pow(t21,2)*t23) + t9);
return se;
}

double specific_momentum(double r)
{
double t1 = pow(chi,2);
double t3 = pow(r,2);
double t4 = pow(spin,2);
double t5 = t1*t4;
double t6 = t3 + t5;
double t15 = pow(spin,4);
double t2 = -1 + t1;
double t29 = -t3;
double t30 = t29 + t5;
double t35 = pow(t6,-4);
double t18 = pow(r,5);
double t19 = 2*r*t1;
double t20 = -1 + t1 + t19;
double t21 = t20*t3*t4;
double t22 = pow(chi,4);
double t23 = -1 + r;
double t24 = t22*t23;
double t25 = t1 + t24;
double t26 = t15*t25;
double t27 = t18 + t21 + t26;
double t7 = 1/t6;
double t28 = 1/t27;
double t31 = pow(t6,-2);
double t32 = -2*spin*t2*t30*t31;
double t33 = pow(t2,2);
double t34 = pow(t30,2);
double t36 = 4*t33*t34*t35*t4;
double t37 = 4*t2*t27*t30*t35;
double t38 = t36 + t37;
double t39 = sqrt(t38);
double t40 = t32 + t39;
double t9 = pow(r,4);
double t11 = -2*t1;
double t12 = r*t1;
double t13 = 2 + r + t11 + t12;
double t14 = r*t13*t4;
double t16 = t1*t15;
double t17 = t14 + t16 + t9;
double Lz = (-0.5*(t17*t28*t40*t6) - 2*r*spin*t2*t7)/sqrt(2*r*spin*t28*t40*t6 + (t17*pow(t40,2)*pow(t6,3))/(4.*t2*pow(t27,2)) + (-2*r + t3 + t5)*t7);
return Lz;
}