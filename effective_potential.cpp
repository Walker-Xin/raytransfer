long double Veff_deri2(long double r, long double E, long double Lz)
{
long double t7 = exp(2);
long double t11 = pow(r,2);
long double t8 = pow(r,9);
long double t9 = -8*r;
long double t12 = 3*t11;
long double t13 = 6 + t12 + t9;
long double t30 = pow(r,7);
long double t15 = pow(defpar,2);
long double t26 = pow(r,4);
long double t21 = pow(Lz,2);
long double t4 = pow(spin,2);
long double t2 = -2 + r;
long double t33 = 15*t11;
long double t57 = pow(r,5);
long double t18 = 21*t11;
long double t58 = pow(r,3);
long double t16 = pow(r,6);
long double d2Veff = (2*(-2*E*Lz*spin*t30*(4*t11*t13 + defpar*(40 - 48*r + t33)) + 2*E*Lz*pow(spin,5)*t26*(5*defpar*(14 - 9*r) - 2*t58) - 30*defpar*E*Lz*pow(spin,7)*t58 - 6*E*Lz*pow(spin,3)*t57*(defpar*(28 - 42*r + t33) + 2*t2*t58) + t15*t16*(60 - 70*r + t18)*t7 + 6*defpar*pow(spin,8)*(6*defpar + 5*t58)*t7 + r*pow(spin,6)*((-160 + 129*r)*t15 + 2*t16 - 2*defpar*(70 - 51*r)*t58)*t7 + t4*((224 - 300*r + 99*t11)*t15*t26*t7 + 2*defpar*(80 - 102*r + 33*t11)*t30*t7 + pow(r,8)*((-12 + 6*r + t11)*t21 + 2*r*(12 - 14*r + t12)*t7)) + pow(spin,4)*(3*t11*(60 - 130*r + 57*t11)*t15*t7 + 6*defpar*(28 - 52*r + t18)*t57*t7 + 2*t30*(t21 + 3*r*t2*t7)) + 4*defpar*t13*t7*t8 + (3*pow(2 - r,3)*t21 + 2*t26*t7)*t8))/(pow(r,10)*pow(r*t2 + t4,3));
return d2Veff;
}