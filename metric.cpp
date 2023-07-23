/* Metric and derivative of metric used throughout code */
/* If possible, use Maple or Mathematica to optimize code */
// TODO: use python to optimise code

void metric(long double r, long double th, long double g[][4])
{
	long double gtt, grr, gthth, gpp, gtp, gtr, grp;

	long double t1 = pow(r,2);
	long double t2 = pow(spin,2);
	long double t3 = cos(th);
	long double t4 = pow(t3,2);
	long double t5 = t2*t4;
	long double t6 = t1 + t5;
	long double t7 = 2*th;
	long double t8 = cos(t7);
	long double t9 = t2*t8;
	long double t11 = t1 + t9;
	long double t12 = pow(t11,-2);
	long double t15 = sin(th);
	long double t16 = pow(t15,2);
	long double t13 = -2 + r;
	long double t14 = r*t13;

	gtt = t12*(t14 + t2 - 4*t16*t2)*t6;
	grr = t6/(-2*r + t1 + t2);
	gthth = t6;
	gpp = t12*t16*(pow(t1 + t2,2) - t16*t2*(t14 + t2))*t6;
	gtp = -(spin*t12*t16*(r*(2 + r) + t2)*t6);
	gtr = (b1*spin)/r;
	grp = 0;

	g[0][0] = gtt;
	g[0][1] = gtr;
	g[0][3] = gtp;
	g[1][0] = g[0][1];
	g[1][1] = grr;
	g[1][3] = grp;
	g[2][2] = gthth;
	g[3][0] = g[0][3];
	g[3][1] = g[1][3];
	g[3][3] = gpp;
}

void metric_rderivatives(long double r, long double th, long double dg[][4])
{
	long double dgttdr, dgrrdr, dgththdr, dgppdr, dgtpdr, dgtrdr, dgrpdr;

	long double t1 = pow(r,2);
	long double t2 = pow(spin,2);
	long double t14 = pow(spin,4);
	long double t3 = 2*th;
	long double t4 = cos(t3);
	long double t24 = -r;
	long double t25 = 1 + t24;
	long double t5 = t2*t4;
	long double t6 = t1 + t5;
	long double t7 = pow(t6,-3);
	long double t8 = pow(r,4);
	long double t12 = pow(r,3);
	long double t26 = 4*th;
	long double t27 = cos(t26);
	long double t64 = pow(spin,6);
	long double t65 = 6*th;
	long double t66 = cos(t65);
	long double t77 = 1 + r;
	long double t70 = sin(th);
	long double t71 = pow(t70,2);

	dgttdr = (t7*(-t14 + 5*r*t14 + 6*t1*t2 + 2*t12*t2 - t14*t25*t27 - 2*t2*((3 + r)*t1 + (1 + 3*r)*t2)*t4 + 4*t8))/2.;
	dgrrdr = (2*(r*(t2 + t24) + t2*t25*pow(cos(th),2)))/pow((-2 + r)*r + t2,2);
	dgththdr = 2*r;
	dgppdr = (t7*t71*(16*pow(r,7) - 18*t1*t14 - 6*t12*t14 + 2*r*t14*(-3*r + 7*t1 + 4*t2)*t27 - t64*t66 + r*t64*t66 - 8*t2*t8 + t2*t4*((1 + 7*r)*t14 + 8*(3 + 5*r)*t1*t2 + 8*(1 + 6*r)*t8)))/8.;
	dgtpdr = -0.5*(spin*t7*t71*(t14 - 3*r*t14 - 6*t1*t2 - 6*t12*t2 + t14*t27*t77 + 2*t2*(3*t1 + t2)*t4*t77 - 4*t8));
	dgtrdr = -((b1*spin)/pow(r,2));
	dgrpdr = 0;

 	dg[0][0] = dgttdr;
	dg[0][1] = dgtrdr;
	dg[0][3] = dgtpdr;
	dg[1][0] = dg[0][1];
	dg[1][1] = dgrrdr;
	dg[1][3] = dgrpdr;
	dg[2][2] = dgththdr;
	dg[3][0] = dg[0][3];
	dg[3][1] = dg[1][3];
	dg[3][3] = dgppdr;
}

void metric_r2derivatives(long double r, long double th, long double dg2[][4])
{
	long double dgttdr2, dgrrdr2, dgththdr2, dgppdr2, dgtpdr2, dgtrdr2, dgrpdr2;

	long double t2 = pow(spin,2);
	long double t1 = pow(r,2);
	long double t11 = pow(r,3);
	long double t13 = pow(r,4);
	long double t15 = pow(spin,4);
	long double t3 = 2*th;
	long double t4 = cos(t3);
	long double t17 = pow(spin,6);
	long double t5 = t2*t4;
	long double t6 = t1 + t5;
	long double t7 = pow(t6,-4);
	long double t8 = pow(r,5);
	long double t29 = 4*th;
	long double t30 = cos(t29);
	long double t32 = 6*th;
	long double t33 = cos(t32);
	long double t57 = pow(spin,8);
	long double t9 = 16*t8;
	long double t12 = 48*t11*t2;
	long double t19 = 80*t11;
	long double t21 = 48*r*t2;
	long double t81 = sin(th);
	long double t82 = pow(t81,2);

	dgttdr2 = -0.25*(t7*(t12 + 56*t1*t15 + 6*t17 + 12*t13*t2 + 2*t15*(8*t1 + 3*t2)*t30 - t17*t33 - t2*(12*t13 + 11*t15 + t19 + 72*t1*t2 + t21)*t4 + t9));
	dgrrdr2 = (2*(2*t11 + t15 - 3*t1*t2 - t2*(-4 + 6*r - 3*t1 + t2)*pow(cos(th),2)))/pow(-2*r + t1 + t2,3);
	dgththdr2 = 2;
	dgppdr2 = (t7*t82*(32*pow(r,8) + 176*t11*t15 + 276*t13*t15 + 48*r*t17 + 120*t1*t17 + 4*t15*(20*t11 + 39*t13 + 2*t15 + 12*r*t2 + 10*t1*t2)*t30 + 32*t1*t17*t33 + 7*t57 + 8*t33*t57 + 8*t2*t4*(16*pow(r,6) - 12*r*t15 - 8*t1*t15 + t17 - 32*t11*t2 - 30*t13*t2 - 4*t8) + 32*t2*t8 + t57*cos(8*th)))/16.;
	dgtpdr2 = -0.25*(spin*t7*t82*(t12 + 48*t1*t15 + 2*t17 + 36*t13*t2 + 2*t15*(4*t1 + t2)*t30 + t17*t33 - t2*(36*t13 + 5*t15 + t19 + 56*t1*t2 + t21)*t4 + t9));
	dgtrdr2 = (2*b1*spin)/pow(r,3);
	dgrpdr2 = 0;

	dg2[0][0] = dgttdr2;
	dg2[0][1] = dgtrdr2;
	dg2[0][3] = dgtpdr2;
	dg2[1][0] = dg2[0][1];
	dg2[1][1] = dgrrdr2;
	dg2[1][3] = dgrpdr2;
	dg2[2][2] = dgththdr2;
	dg2[3][0] = dg2[0][3];
	dg2[3][1] = dg2[1][3];
	dg2[3][3] = dgppdr2;
}

void uppermetric(long double r, long double th, long double rth[])
{
	long double gutt, gurr, guthth, gupp, gutp, gurp;

	long double t1 = pow(r,2);
	long double t2 = pow(spin,2);
	long double t14 = sin(th);
	long double t15 = pow(t14,2);
	long double t3 = cos(th);
	long double t4 = pow(t3,2);
	long double t5 = t2*t4;
	long double t6 = t1 + t5;
	long double t16 = t1 + t2;
	long double t17 = pow(t16,2);
	long double t18 = -2 + r;
	long double t19 = r*t18;
	long double t20 = t19 + t2;
	long double t21 = -(t15*t2*t20);
	long double t22 = t17 + t21;
	long double t8 = 2*th;
	long double t9 = cos(t8);
	long double t11 = t2*t9;
	long double t12 = t1 + t11;
	long double t13 = pow(t12,2);
	long double t7 = pow(t6,2);
	long double t23 = 2 + r;
	long double t24 = r*t23;
	long double t25 = t2 + t24;
	long double t26 = pow(t25,2);
	long double t36 = -4*t15*t2;
	long double t37 = t19 + t2 + t36;
	long double t27 = pow(t6,3);
	long double t28 = pow(t14,4);
	long double t29 = t1*t2*t26*t27*t28;
	long double t30 = pow(b1,2);
	long double t31 = t13*t2*t20*t30;
	long double t32 = pow(r,3);
	long double t33 = r*t2*t4;
	long double t34 = t32 + t33;
	long double t35 = pow(t34,2);
	long double t38 = -(t35*t37);
	long double t39 = t31 + t38;
	long double t40 = t15*t22*t39*t6;
	long double t41 = t29 + t40;
	long double t42 = 1/t41;
	
	gutt = -(t1*t13*t15*t22*t42*t7);
	gurr = -(t1*t15*t20*(-(t15*t2*t26) + t22*t37)*t42*t7);
	guthth = 1/t6;
	gupp = -(t1*pow(t12,4)*t20*t42*(-((t2*t30)/pow(r,2)) + (t37*t7)/(pow(t12,2)*t20)));
	gutp = -(spin*t1*t13*t15*t25*t42*t7);
	gurp = b1*r*t13*t15*t2*t20*t25*t42*t6;

	rth[0] = gurr;
	rth[1] = guthth;
}
