double Veff_deri2(double r, double E, double Lz)
{
double t7 = pow(r,2);
double t8 = pow(Lz,2);
double t15 = exp(2);
double t21 = 3*t7;
double t4 = pow(spin,2);
double t9 = -r;
double t12 = 2 + t9;
double t2 = -2 + r;
double d2Veff = (2*(-4*E*Lz*pow(spin,5) + 12*E*Lz*r*pow(spin,3)*t12 + 2*pow(spin,6)*t15 - 8*E*Lz*spin*(6 - 8*r + t21)*t7 + 2*pow(spin,4)*(3*r*t15*t2 + t8) + t7*(2*pow(r,4)*t15 + 3*pow(t12,3)*t8) + r*t4*(2*r*t15*(12 - 14*r + t21) + (-12 + 6*r + t7)*t8)))/(pow(r,3)*pow(r*t2 + t4,3));
return d2Veff;
}