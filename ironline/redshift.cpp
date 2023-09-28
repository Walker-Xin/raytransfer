void redshift(double r, double th, double ktkp, double &gg)
{
double t3 = pow(spin,2);
double t1 = pow(r,2);
double t15 = pow(spin,4);
double t4 = 2*th;
double t5 = cos(t4);
double t6 = t3*t5;
double t2 = 2*t1;
double t7 = t2 + t3 + t6;
double t31 = -2*t1;
double t32 = t3 + t31 + t6;
double t9 = pow(r,5);
double t11 = 8*t9;
double t12 = -4*t1*t3;
double t13 = pow(r,3);
double t14 = 8*t13*t3;
double t16 = 3*r*t15;
double t17 = 2*r;
double t18 = 1 + t17;
double t19 = r*t18;
double t20 = t19 + t3;
double t21 = 4*r*t20*t3*t5;
double t22 = -r;
double t23 = 1 + t22;
double t24 = 4*th;
double t25 = cos(t24);
double t26 = -(t15*t23*t25);
double t27 = t11 + t12 + t14 + t15 + t16 + t21 + t26;
double t34 = sin(th);
double t35 = pow(t34,2);
double t37 = pow(t7,-4);
double t28 = 1/t27;
double t33 = pow(t7,-2);
double t36 = 4*spin*t32*t33*t35;
double t38 = -4*t27*t32*t35*t37;
double t39 = pow(t32,2);
double t40 = pow(t34,4);
double t41 = 16*t3*t37*t39*t40;
double t42 = t38 + t41;
double t43 = sqrt(t42);
double t44 = t36 + t43;
double t29 = 1/sin(th);
double t30 = pow(t29,2);
double t48 = -2 + r;
gg = sqrt((t3 + 2*r*t48 + t6)/t7 + 8*r*spin*t28*t44*t7 - (t30*pow(t44,2)*(pow(t1 + t3,2) - t3*t35*(t3 + r*t48))*pow(t7,4)*(t1 + t3*pow(cos(th),2)))/(pow(t27,2)*pow(t1 + t3 - t3*t35,2)))/(1 - ktkp*t28*t30*t44*pow(t7,2));
}

void redshift_bambi(double spin, double spin2, double epsilon_r, double epsilon_t, double radius, double ktt, double ktkp, double kyy, double& gg, double& ldr)
{
	double r, r2, r3;
	double dr, temp;
	double cc;
	double def_r,def_t;
	double g00[3], g03[3], g33[3];
	double g001, g031, g331;
	double Omega;
	double uet;
	double uephi;
	double mem;
	double H;
	
	int i;
	
	
	dr   = 0.001*radius;
	temp = radius + dr;
	dr   = temp - radius;
	
	cc = ktkp;
	
	for (i = 0; i <= 2; i++) {
		
		r  = radius + dr*(i - 1);
		r2 = r*r;
		r3 = r2*r;
		
		def_t = epsilon_t/r3;
		def_r = epsilon_r/r3;
		H = (1+def_r)*(1+def_t);
		H = sqrt(H);

		g00[i] = - (1 + def_t)*(1 - 2/r);
		g03[i] = - spin*(H-(1 + def_t)*(1 - 2/r));
		g33[i] = r2 + spin2*(2*H- (1 + def_t)*(1 - 2/r));
		
	}
	
	g001 = 0.5*(g00[2] - g00[0])/dr;
	g031 = 0.5*(g03[2] - g03[0])/dr;
	g331 = 0.5*(g33[2] - g33[0])/dr;
	
	Omega  = (-g031 + sqrt(g031*g031 - g001*g331))/g331;
	
	uet = sqrt(-g00[1] - 2*g03[1]*Omega - g33[1]*Omega*Omega);
	
	gg = uet/(1 - cc*Omega);
	
	uet = 1/uet;
	
	uephi = Omega*uet;
	
	mem = kyy*radius;
	mem = mem/(ktt*uet + ktt*cc*uephi);
	
	if (mem < 0) mem = 0;
	if (mem > 1) mem = 1;
	
	ldr = 0.5 + 0.75*mem;

	return ;
}