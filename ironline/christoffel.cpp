void christoffel(double r, double th, double christ[4][4][4])
{
double t3 = pow(spin,2);
double t1 = pow(r,2);
double t4 = 2*th;
double t5 = cos(t4);
double t6 = t3*t5;
double t12 = pow(r,4);
double t26 = pow(defpar,2);
double t17 = pow(spin,4);
double t18 = -2 + r;
double t19 = r*t18;
double t20 = t19 + t3;
double t8 = 2*t1;
double t23 = pow(r,6);
double t24 = 8*t23;
double t25 = 8*t12*t3;
double t27 = 2*t12*t26*t3;
double t28 = 2*r*t17*t26;
double t29 = 3*t1*t17;
double t30 = 3*t1*t17*t26;
double t31 = pow(spin,6);
double t32 = t26*t31;
double t33 = 8*t12;
double t34 = t18*t26;
double t35 = 4*r;
double t36 = t34 + t35;
double t37 = r*t3*t36;
double t38 = t17*t26;
double t39 = t33 + t37 + t38;
double t40 = t3*t39*t5;
double t41 = 4*th;
double t42 = cos(t41);
double t43 = t1*t17*t42;
double t44 = t24 + t25 + t27 + t28 + t29 + t30 + t32 + t40 + t43;
double t45 = 1/t44;
double t2 = -2*t1;
double t7 = t2 + t3 + t6;
double t9 = t3 + t6 + t8;
double t11 = pow(t9,-2);
double t13 = 2*t12;
double t14 = 3*r;
double t15 = 2 + t14;
double t16 = r*t15*t3;
double t21 = t20*t3*t5;
double t22 = t13 + t16 + t17 + t21;
double t47 = 1/t20;
double t48 = t1 + t3;
double t49 = -t3;
double t50 = -(t3*t5);
double t51 = t49 + t50 + t8;
double t52 = 4*t1*t45*t47*t48*t51;
double t53 = pow(r,3);
double t54 = sin(t4);
double t56 = sin(th);
double t57 = pow(t56,2);
double t55 = -8*t3*t45*t53*t54;
double t93 = pow(spin,3);
double t94 = -(defpar*r*t22*t45*t47*t54*t93);
double t58 = 4*defpar*r*t11*t22*t3*t45*t57*t7;
double t95 = -6*t12;
double t96 = -3*t1*t3;
double t97 = -t1;
double t98 = t3 + t97;
double t99 = t3*t5*t98;
double t100 = t17 + t95 + t96 + t99;
double t101 = 4*spin*t1*t100*t45*t47*t57;
double t103 = cos(th);
double t104 = pow(t56,3);
double t105 = 16*t103*t104*t45*t53*t93;
double t65 = pow(r,5);
double t84 = -r;
double t85 = 2 + t84;
double t121 = 1/t9;
double t123 = 4*defpar*r*spin*t121*t45*t48*t7;
double t60 = pow(r,7);
double t73 = pow(spin,8);
double t115 = 1 + t84;
double t124 = 8*defpar*t1*t121*t20*t45*t54*t93;
double t150 = -2*t1*t3*t45*t54*t9;
double t125 = 8*spin*t1*t121*t20*t45*t57*t7;
double t151 = -4*defpar*r*t100*t121*t3*t45*t57;
double t153 = -16*defpar*t1*t103*t104*t121*t17*t20*t45;
double t106 = 8*t65;
double t107 = -4*t1*t3;
double t108 = 8*t3*t53;
double t109 = 3*r*t17;
double t110 = 2*r;
double t111 = 1 + t110;
double t112 = r*t111;
double t113 = t112 + t3;
double t114 = 4*r*t113*t3*t5;
double t116 = -(t115*t17*t42);
double t117 = t106 + t107 + t108 + t109 + t114 + t116 + t17;
double t155 = pow(t9,-3);
double t158 = pow(t103,2);
double t159 = t158*t3;
double t160 = t1 + t159;
double t161 = 1/t160;
double t163 = r*t161;
double t157 = 8*r*spin*t155*t48*t54;
double t178 = 4*t12;
double t181 = 2*t1*t3*t5;
double t179 = t110 + t34;
double t180 = r*t179*t3;
double t182 = t178 + t180 + t181 + t38;
double t183 = -2*spin*t121*t182*t45*t47*t7;
double t59 = pow(t20,-2);
double t87 = -3 + r;
double t80 = -1 + t35;
double t184 = 2 + t26;
double t185 = t1*t184*t3;
double t186 = t178 + t181 + t185 + t38;
double t187 = 1/tan(th);
double t188 = -8*r*spin*t121*t186*t187*t45;
double t196 = -8*defpar*t1*t103*t17*t45*t47*t56;
double t62 = pow(r,8);
double t219 = -2 + t14;
double t221 = 6*r;
double t222 = t221 + t34;
double t189 = 16*defpar*t1*t11*t45*t57*t7*t93;
double t197 = 1/sin(th);
double t198 = pow(t197,2);
double t199 = 4*t18*t53;
double t200 = t180 + t181 + t199 + t38;
double t201 = t117*t200*t57*t9;
double t202 = pow(t56,4);
double t203 = 64*t160*t202*t3*t53*t7;
double t204 = t201 + t203;
double t205 = (t11*t198*t204*t45*t47)/2.;
double t207 = 32*t62;
double t208 = 32*t3*t65;
double t209 = 48*t23*t3;
double t210 = 8*t23*t26*t3;
double t211 = 8*t17*t53;
double t212 = 16*t17*t26*t53;
double t213 = 36*t12*t17;
double t214 = 16*t12*t17*t26;
double t215 = 10*r*t26*t31;
double t216 = 10*t1*t31;
double t217 = 11*t1*t26*t31;
double t218 = 3*t26*t73;
double t220 = 16*t219*t65;
double t223 = 8*t222*t3*t53;
double t224 = 15*r;
double t225 = 12*r;
double t226 = -8 + t225;
double t227 = t226*t26;
double t228 = t224 + t227;
double t229 = r*t17*t228;
double t230 = 4*t26*t31;
double t231 = t220 + t223 + t229 + t230;
double t232 = t231*t3*t5;
double t233 = 4*t219*t53;
double t234 = r*t222*t3;
double t235 = t233 + t234 + t38;
double t236 = t17*t235*t42;
double t237 = 6*th;
double t238 = cos(t237);
double t239 = t1*t238*t31;
double t240 = t207 + t208 + t209 + t210 + t211 + t212 + t213 + t214 + t215 + t216 + t217 + t218 + t232 + t236 + t239;
double t241 = (t121*t187*t240*t45)/2.;
christ[0][0][0] = -4*defpar*r*spin*t11*t22*t45*t7;
christ[0][0][1] = t52;
christ[0][0][2] = t55;
christ[0][0][3] = t58;
christ[0][1][0] = t52;
christ[0][1][1] = (defpar*spin*t45*t59*(-2*t1*t17 + 32*t12*t17 - 24*t12*t3 + 28*t23*t3 - r*t31 + 15*t1*t31 - 3*t17*t53 - 24*t60 + 8*t62 - 32*t3*t65 + 3*t73 + 4*t3*t5*((6 - 6*r + t1)*t12 + t31 + t3*(-9 + t35)*t53 + r*t17*t80) + t17*t42*(t17 + t1*t85 + r*t3*t87)))/2.;
christ[0][1][2] = t94;
christ[0][1][3] = t101;
christ[0][2][0] = t55;
christ[0][2][1] = t94;
christ[0][2][2] = -2*defpar*spin*t1*t22*t45;
christ[0][2][3] = t105;
christ[0][3][0] = t58;
christ[0][3][1] = t101;
christ[0][3][2] = t105;
christ[0][3][3] = -(defpar*r*spin*t11*t117*t22*t45*t57);
christ[1][0][0] = -8*t1*t121*t45*t51*(t49 + r*t85);
christ[1][0][1] = t123;
christ[1][0][2] = t124;
christ[1][0][3] = t125;
christ[1][1][0] = t123;
christ[1][1][1] = (t45*t47*(t12*t17 + 4*t1*t17*t26 - 5*t12*t17*t26 + 4*t23*t3 - 2*t23*t26*t3 - 4*t1*t26*t31 - t3*(4*t23 + 2*r*t17*t18*t26 + t1*(-4*r + pow(t18,2)*t26)*t3 + t32)*t5 + 3*t17*t53 + 4*t17*t26*t53 + t115*t17*t42*t53 - 8*t60 + 4*t26*t3*t65 - t26*t73))/r;
christ[1][1][2] = t150;
christ[1][1][3] = t151;
christ[1][2][0] = t124;
christ[1][2][1] = t150;
christ[1][2][2] = -4*t20*t45*t53*t9;
christ[1][2][3] = t153;
christ[1][3][0] = t125;
christ[1][3][1] = t151;
christ[1][3][2] = t153;
christ[1][3][3] = -2*t1*t117*t121*t20*t45*t57;
christ[2][0][0] = -8*r*t155*t3*t54;
christ[2][0][1] = 0;
christ[2][0][2] = 0;
christ[2][0][3] = t157;
christ[2][1][0] = 0;
christ[2][1][1] = t103*t161*t3*t47*t56;
christ[2][1][2] = t163;
christ[2][1][3] = 0;
christ[2][2][0] = 0;
christ[2][2][1] = t163;
christ[2][2][2] = -(t103*t161*t3*t56);
christ[2][2][3] = 0;
christ[2][3][0] = t157;
christ[2][3][1] = 0;
christ[2][3][2] = 0;
christ[2][3][3] = -0.5*(t155*t54*(10*r*t17 + 11*t1*t17 + t24 + 16*t12*t3 + 3*t31 + t17*(-2*r + t1 + t3)*t42 + 16*t3*t53 + 4*t20*t3*t5*(t3 + t8)));
christ[3][0][0] = -16*defpar*t1*t11*t3*t45*t7;
christ[3][0][1] = t183;
christ[3][0][2] = t188;
christ[3][0][3] = t189;
christ[3][1][0] = t183;
christ[3][1][1] = 4*defpar*r*t3*t45*t59*(t17 + r*t3*t80 + t3*t5*(t3 + t84) + 2*t53*t87);
christ[3][1][2] = t196;
christ[3][1][3] = t205;
christ[3][2][0] = t188;
christ[3][2][1] = t196;
christ[3][2][2] = -8*defpar*t3*t45*t53;
christ[3][2][3] = t241;
christ[3][3][0] = t189;
christ[3][3][1] = t205;
christ[3][3][2] = t241;
christ[3][3][3] = 4*defpar*t1*t11*t3*t45*t57*(-t17 - 3*r*t17 + 4*t1*t3 + t115*t17*t42 - 4*r*t113*t3*t5 - 8*t3*t53 - 8*t65);
}

void Christoffel_jiale(double w1, double w2, double CS[][4][4])
{

	double g[4][4];
	double gg;
	double gr[4][4];
	double grr[4][4];
	double gl[4][4];
	double gll[4][4];
	double dw1, dw2, temp;
	
	double invg[4][4];
	double Dg[4][4][4];
	
	
	/* ----- set dw1 and dw2 ----- */
	
	dw1  = 0.01*w1;
	temp = w1 + dw1;
	dw1  = temp - w1;
	
	dw2  = 0.01;
	temp = w2 + dw2;
	dw2  = temp - w2;
	
	
	/* ----- compute inverse metric ----- */
	
	metric(w1, w2, g);
	

/*	gg = g[0][0]*g[3][3] - g[0][3]*g[0][3];
	
	invg[0][0] = g[3][3]/gg; 
	invg[0][3] = - g[0][3]/gg;
	invg[1][1] = 1/g[1][1];
	invg[2][2] = 1/g[2][2];
	invg[3][0] = invg[0][3]; 
	invg[3][3] = g[0][0]/gg;*/

	gg = (g[0][1]*g[0][1] - g[0][0]*g[1][1])*g[3][3]+g[0][3]*g[0][3]*g[1][1];
	
	invg[0][0] = -g[1][1]*g[3][3]/gg; 
	invg[0][1] = g[0][1]*g[3][3]/gg;
	invg[0][3] = g[0][3]*g[1][1]/gg;
	invg[1][0] = invg[0][1];
	invg[1][1] = (g[0][3]*g[0][3]-g[0][0]*g[3][3])/gg;
	invg[1][3] = -g[0][1]*g[0][3]/gg;
	invg[2][2] = 1/g[2][2];
	invg[3][0] = invg[0][3];
	invg[3][1] = invg[1][3];
	invg[3][3] = (g[0][1]*g[0][1] - g[0][0]*g[1][1])/gg;

/*
	gg = (g[1][3]*g[3][1] - g[1][1]*g[3][3]) * g[0][0] + g[0][3]*g[3][0]*g[1][1];
	
	invg[0][0] = (g[1][3]*g[3][1] - g[1][1]*g[3][3])/gg;
	invg[0][1] = -g[0][3]*g[3][1]/gg;
	invg[0][3] = g[0][3]*g[1][1]/gg;
	invg[1][0] = invg[0][1];
	invg[1][1] = (g[0][3]*g[3][0] - g[0][0]*g[3][3])/gg;
	invg[1][3] = g[0][0]*g[1][3]/gg;
	invg[2][2] = 1/g[2][2];
	invg[3][0] = invg[0][3];
	invg[3][1] = invg[1][3];
	invg[3][3] = -g[0][0]*g[1][1]/gg;
*/
	
		
	/* ----- compute metric derivatives ----- */
	
	metric(w1+dw1,w2,gr);
	metric(w1-dw1,w2,gl);
	
	Dg[0][0][1] = 0.5*(gr[0][0] - gl[0][0])/dw1;
	Dg[0][3][1] = 0.5*(gr[0][3] - gl[0][3])/dw1;
	Dg[1][1][1] = 0.5*(gr[1][1] - gl[1][1])/dw1;
	Dg[2][2][1] = 0.5*(gr[2][2] - gl[2][2])/dw1;
	Dg[3][0][1] = Dg[0][3][1];
	Dg[3][3][1] = 0.5*(gr[3][3] - gl[3][3])/dw1;

	 	
	metric(w1,w2+dw2,gr);
	metric(w1,w2-dw2,gl);
	Dg[0][0][2] = 0.5*(gr[0][0] - gl[0][0])/dw2;
	Dg[0][3][2] = 0.5*(gr[0][3] - gl[0][3])/dw2;
	Dg[1][1][2] = 0.5*(gr[1][1] - gl[1][1])/dw2;
	Dg[2][2][2] = 0.5*(gr[2][2] - gl[2][2])/dw2;
	Dg[3][0][2] = Dg[0][3][2];
	Dg[3][3][2] = 0.5*(gr[3][3] - gl[3][3])/dw2;
	
		
	
	/*
	metric(w1+dw1,w2,gr);
	metric(w1+dw1+dw1,w2,grr);
	metric(w1-dw1,w2,gl);
	metric(w1-dw1-dw1,w2,gll);
	
	Dg[0][0][1] = (8*gr[0][0] - grr[0][0] - 8*gl[0][0] + gll[0][0])/dw1/12;
	Dg[0][3][1] = (8*gr[0][3] - grr[0][3] - 8*gl[0][3] + gll[0][3])/dw1/12;
	Dg[1][1][1] = (8*gr[1][1] - grr[1][1] - 8*gl[1][1] + gll[1][1])/dw1/12;
	Dg[2][2][1] = (8*gr[2][2] - grr[2][2] - 8*gl[2][2] + gll[2][2])/dw1/12;
	Dg[3][0][1] = Dg[0][3][1];
	Dg[3][3][1] = (8*gr[3][3] - grr[3][3] - 8*gl[3][3] + gll[3][3])/dw1/12;
	
	metric(w1,w2+dw2,gr);
	metric(w1,w2+dw2+dw2,grr);
	metric(w1,w2-dw2,gl);
	metric(w1,w2-dw2-dw2,gll);
	
	Dg[0][0][2] = (8*gr[0][0] - grr[0][0] - 8*gl[0][0] + gll[0][0])/dw2/12;
	Dg[0][3][2] = (8*gr[0][3] - grr[0][3] - 8*gl[0][3] + gll[0][3])/dw2/12;
	Dg[1][1][2] = (8*gr[1][1] - grr[1][1] - 8*gl[1][1] + gll[1][1])/dw2/12;
	Dg[2][2][2] = (8*gr[2][2] - grr[2][2] - 8*gl[2][2] + gll[2][2])/dw2/12;
	Dg[3][0][2] = Dg[0][3][2];
	Dg[3][3][2] = (8*gr[3][3] - grr[3][3] - 8*gl[3][3] + gll[3][3])/dw2/12;
	*/
	
	/* ----- compute Christoffel symbols ----- */
	
/*	CS[0][0][0] = 0;
	CS[0][0][1] = invg[0][0]*Dg[0][0][1] + invg[0][3]*Dg[0][3][1];
	CS[0][0][2] = invg[0][0]*Dg[0][0][2] + invg[0][3]*Dg[0][3][2];
	CS[0][0][3] = 0;
	CS[0][1][0] = CS[0][0][1];
	CS[0][1][1] = 0;
	CS[0][1][2] = 0;
	CS[0][1][3] = invg[0][0]*Dg[0][3][1] + invg[0][3]*Dg[3][3][1];
	CS[0][2][0] = CS[0][0][2];
	CS[0][2][1] = 0;
	CS[0][2][2] = 0;
	CS[0][2][3] = invg[0][0]*Dg[0][3][2] + invg[0][3]*Dg[3][3][2];
	CS[0][3][0] = 0;
	CS[0][3][1] = CS[0][1][3];
	CS[0][3][2] = CS[0][2][3];
	CS[0][3][3] = 0;
	
	CS[1][0][0] = - invg[1][1]*Dg[0][0][1];
	CS[1][0][1] = 0;
	CS[1][0][2] = 0;
	CS[1][0][3] = - invg[1][1]*Dg[0][3][1];
	CS[1][1][0] = 0;
	CS[1][1][1] = invg[1][1]*Dg[1][1][1];
	CS[1][1][2] = invg[1][1]*Dg[1][1][2];
	CS[1][1][3] = 0;
	CS[1][2][0] = 0;
	CS[1][2][1] = CS[1][1][2];
	CS[1][2][2] = - invg[1][1]*Dg[2][2][1];
	CS[1][2][3] = 0;
	CS[1][3][0] = CS[1][0][3];
	CS[1][3][1] = 0;
	CS[1][3][2] = 0;
	CS[1][3][3] = - invg[1][1]*Dg[3][3][1];
	
	CS[2][0][0] = - invg[2][2]*Dg[0][0][2];
	CS[2][0][1] = 0;
	CS[2][0][2] = 0;
	CS[2][0][3] = - invg[2][2]*Dg[0][3][2];
	CS[2][1][0] = 0;
	CS[2][1][1] = - invg[2][2]*Dg[1][1][2];
	CS[2][1][2] = invg[2][2]*Dg[2][2][1];
	CS[2][1][3] = 0;
	CS[2][2][0] = 0;
	CS[2][2][1] = CS[2][1][2];
	CS[2][2][2] = invg[2][2]*Dg[2][2][2];
	CS[2][2][3] = 0;
	CS[2][3][0] = CS[2][0][3];
	CS[2][3][1] = 0;
	CS[2][3][2] = 0;
	CS[2][3][3] = - invg[2][2]*Dg[3][3][2];
	
	CS[3][0][0] = 0;
	CS[3][0][1] = invg[3][3]*Dg[0][3][1] + invg[3][0]*Dg[0][0][1];
	CS[3][0][2] = invg[3][3]*Dg[0][3][2] + invg[3][0]*Dg[0][0][2];
	CS[3][0][3] = 0;
	CS[3][1][0] = CS[3][0][1];
	CS[3][1][1] = 0;
	CS[3][1][2] = 0;
	CS[3][1][3] = invg[3][3]*Dg[3][3][1] + invg[3][0]*Dg[0][3][1];
	CS[3][2][0] = CS[3][0][2];
	CS[3][2][1] = 0;
	CS[3][2][2] = 0;
	CS[3][2][3] = invg[3][3]*Dg[3][3][2] + invg[3][0]*Dg[0][3][2];
	CS[3][3][0] = 0;
	CS[3][3][1] = CS[3][1][3];
	CS[3][3][2] = CS[3][2][3];
	CS[3][3][3] = 0;*/

/*	CS[0][0][0] = -invg[0][1]*Dg[0][0][1];
	CS[0][0][1] = invg[0][0]*Dg[0][0][1] + invg[0][3]*Dg[0][3][1];
	CS[0][0][2] = invg[0][0]*Dg[0][0][2] + invg[0][3]*Dg[0][3][2];
	CS[0][0][3] = -invg[0][1]*Dg[0][3][1];
	CS[0][1][0] = CS[0][0][1];
	CS[0][1][1] = 2*invg[0][3]*Dg[1][3][1] + invg[0][1]*Dg[1][1][1];
	CS[0][1][2] = invg[0][1]*Dg[1][1][2] + invg[0][3]*Dg[1][3][2];
	CS[0][1][3] = invg[0][0]*Dg[0][3][1] + invg[0][3]*Dg[3][3][1];
	CS[0][2][0] = CS[0][0][2];
	CS[0][2][1] = CS[0][1][2];
	CS[0][2][2] = -invg[0][1]*Dg[2][2][1];
	CS[0][2][3] = invg[0][0]*Dg[0][3][2] + invg[0][3]*Dg[3][3][2] + invg[0][1]*Dg[1][3][2];
	CS[0][3][0] = CS[0][0][3];
	CS[0][3][1] = CS[0][1][3];
	CS[0][3][2] = CS[0][2][3];
	CS[0][3][3] = -invg[0][1]*Dg[3][3][1];
	
	CS[1][0][0] = - invg[1][1]*Dg[0][0][1];
	CS[1][0][1] = invg[1][0]*Dg[0][0][1] + invg[1][3]*Dg[0][3][1];
	CS[1][0][2] = invg[1][0]*Dg[0][0][2] + invg[1][3]*Dg[0][3][2];
	CS[1][0][3] = - invg[1][1]*Dg[0][3][1];
	CS[1][1][0] = CS[1][0][1];
	CS[1][1][1] = 2*invg[1][3]*Dg[1][3][1] + invg[1][1]*Dg[1][1][1];
	CS[1][1][2] = invg[1][1]*Dg[1][1][2] + invg[1][3]*Dg[1][3][2];
	CS[1][1][3] = invg[1][0]*Dg[0][3][1] + invg[1][3]*Dg[3][3][1];
	CS[1][2][0] = CS[1][0][2];
	CS[1][2][1] = CS[1][1][2];
	CS[1][2][2] = - invg[1][1]*Dg[2][2][1];
	CS[1][2][3] = invg[1][0]*Dg[0][3][2] + invg[1][3]*Dg[3][3][2] + invg[1][1]*Dg[1][3][2];
	CS[1][3][0] = CS[1][0][3];
	CS[1][3][1] = CS[1][1][3];
	CS[1][3][2] = CS[1][2][3];
	CS[1][3][3] = - invg[1][1]*Dg[3][3][1];
	
	CS[2][0][0] = - invg[2][2]*Dg[0][0][2];
	CS[2][0][1] = 0;
	CS[2][0][2] = 0;
	CS[2][0][3] = - invg[2][2]*Dg[0][3][2];
	CS[2][1][0] = 0;
	CS[2][1][1] = - invg[2][2]*Dg[1][1][2];
	CS[2][1][2] = invg[2][2]*Dg[2][2][1];
	CS[2][1][3] = - invg[2][2]*Dg[1][3][2];
	CS[2][2][0] = 0;
	CS[2][2][1] = CS[2][1][2];
	CS[2][2][2] = invg[2][2]*Dg[2][2][2];
	CS[2][2][3] = 0;
	CS[2][3][0] = CS[2][0][3];
	CS[2][3][1] = CS[2][1][3];
	CS[2][3][2] = 0;
	CS[2][3][3] = - invg[2][2]*Dg[3][3][2];
	
	CS[3][0][0] = -invg[3][1]*Dg[0][0][1];
	CS[3][0][1] = invg[3][3]*Dg[0][3][1] + invg[3][0]*Dg[0][0][1];
	CS[3][0][2] = invg[3][3]*Dg[0][3][2] + invg[3][0]*Dg[0][0][2];
	CS[3][0][3] = -invg[3][1]*Dg[0][3][1];
	CS[3][1][0] = CS[3][0][1];
	CS[3][1][1] = 2*invg[3][3]*Dg[1][3][1] + invg[3][1]*Dg[1][1][1];
	CS[3][1][2] = invg[3][1]*Dg[1][1][2] + invg[3][3]*Dg[1][3][2];
	CS[3][1][3] = invg[3][3]*Dg[3][3][1] + invg[3][0]*Dg[0][3][1];
	CS[3][2][0] = CS[3][0][2];
	CS[3][2][1] = CS[3][1][2];
	CS[3][2][2] = -invg[3][1]*Dg[2][2][1];
	CS[3][2][3] = invg[3][3]*Dg[3][3][2] + invg[3][1]*Dg[1][3][2] + invg[3][0]*Dg[0][3][2];
	CS[3][3][0] = CS[3][0][3];
	CS[3][3][1] = CS[3][1][3];
	CS[3][3][2] = CS[3][2][3];
	CS[3][3][3] = -invg[3][1]*Dg[3][3][1];*/

	CS[0][0][0] = -invg[0][1]*Dg[0][0][1];
	CS[0][0][1] = invg[0][0]*Dg[0][0][1] + invg[0][3]*Dg[0][3][1];
	CS[0][0][2] = invg[0][0]*Dg[0][0][2] + invg[0][3]*Dg[0][3][2];
	CS[0][0][3] = -invg[0][1]*Dg[0][3][1];
	CS[0][1][0] = CS[0][0][1];
	CS[0][1][1] = 2*invg[0][0]*Dg[0][1][1] + invg[0][1]*Dg[1][1][1];
	CS[0][1][2] = invg[0][1]*Dg[1][1][2];
	CS[0][1][3] = invg[0][0]*Dg[0][3][1] + invg[0][3]*Dg[3][3][1];
	CS[0][2][0] = CS[0][0][2];
	CS[0][2][1] = CS[0][1][2];
	CS[0][2][2] = -invg[0][1]*Dg[2][2][1];
	CS[0][2][3] = invg[0][0]*Dg[0][3][2] + invg[0][3]*Dg[3][3][2];
	CS[0][3][0] = CS[0][0][3];
	CS[0][3][1] = CS[0][1][3];
	CS[0][3][2] = CS[0][2][3];
	CS[0][3][3] = -invg[0][1]*Dg[3][3][1];
	
	CS[1][0][0] = - invg[1][1]*Dg[0][0][1];
	CS[1][0][1] = invg[1][0]*Dg[0][0][1] + invg[1][3]*Dg[0][3][1];
	CS[1][0][2] = invg[1][0]*Dg[0][0][2] + invg[1][3]*Dg[0][3][2];
	CS[1][0][3] = - invg[1][1]*Dg[0][3][1];
	CS[1][1][0] = CS[1][0][1];
	CS[1][1][1] = 2*invg[1][0]*Dg[0][1][1] + invg[1][1]*Dg[1][1][1];
	CS[1][1][2] = invg[1][1]*Dg[1][1][2];
	CS[1][1][3] = invg[1][0]*Dg[0][3][1] + invg[1][3]*Dg[3][3][1];
	CS[1][2][0] = CS[1][0][2];
	CS[1][2][1] = CS[1][1][2];
	CS[1][2][2] = - invg[1][1]*Dg[2][2][1];
	CS[1][2][3] = invg[1][0]*Dg[0][3][2] + invg[1][3]*Dg[3][3][2];
	CS[1][3][0] = CS[1][0][3];
	CS[1][3][1] = CS[1][1][3];
	CS[1][3][2] = CS[1][2][3];
	CS[1][3][3] = - invg[1][1]*Dg[3][3][1];
	
	CS[2][0][0] = - invg[2][2]*Dg[0][0][2];
	CS[2][0][1] = 0;
	CS[2][0][2] = 0;
	CS[2][0][3] = - invg[2][2]*Dg[0][3][2];
	CS[2][1][0] = 0;
	CS[2][1][1] = - invg[2][2]*Dg[1][1][2];
	CS[2][1][2] = invg[2][2]*Dg[2][2][1];
	CS[2][1][3] = 0;
	CS[2][2][0] = 0;
	CS[2][2][1] = CS[2][1][2];
	CS[2][2][2] = invg[2][2]*Dg[2][2][2];
	CS[2][2][3] = 0;
	CS[2][3][0] = CS[2][0][3];
	CS[2][3][1] = 0;
	CS[2][3][2] = 0;
	CS[2][3][3] = - invg[2][2]*Dg[3][3][2];
	
	CS[3][0][0] = -invg[3][1]*Dg[0][0][1];
	CS[3][0][1] = invg[3][3]*Dg[0][3][1] + invg[3][0]*Dg[0][0][1];
	CS[3][0][2] = invg[3][3]*Dg[0][3][2] + invg[3][0]*Dg[0][0][2];
	CS[3][0][3] = -invg[0][1]*Dg[0][3][1];
	CS[3][1][0] = CS[3][0][1];
	CS[3][1][1] = 2*invg[3][0]*Dg[0][1][1] + invg[3][1]*Dg[1][1][1];
	CS[3][1][2] = invg[3][1]*Dg[1][1][2];
	CS[3][1][3] = invg[3][3]*Dg[3][3][1] + invg[3][0]*Dg[0][3][1];
	CS[3][2][0] = CS[3][0][2];
	CS[3][2][1] = CS[3][1][2];
	CS[3][2][2] = -invg[3][1]*Dg[2][2][1];
	CS[3][2][3] = invg[3][3]*Dg[3][3][2] + invg[3][0]*Dg[0][3][2];
	CS[3][3][0] = CS[3][0][3];
	CS[3][3][1] = CS[3][1][3];
	CS[3][3][2] = CS[3][2][3];
	CS[3][3][3] = -invg[3][1]*Dg[3][3][1];

	// Divide everything by 2
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			for (int k=0; k<4; k++)
				CS[i][j][k] *= 0.5;

	return ;
}