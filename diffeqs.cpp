void christoffel(long double r, long double th, long double christ[4][4][4])
{
	long double t2 = pow(r,2);
	long double t6 = pow(spin,2);
	long double t25 = pow(spin,4);
	long double t15 = 2*th;
	long double t16 = cos(t15);
	long double t3 = -2 + r;
	long double t5 = r*t3;
	long double t7 = t5 + t6;
	long double t40 = sin(th);
	long double t41 = pow(t40,2);
	long double t11 = cos(th);
	long double t12 = pow(t11,2);
	long double t50 = t12*t6;
	long double t51 = t2 + t50;
	long double t42 = t2 + t6;
	long double t43 = pow(t42,2);
	long double t44 = -(t41*t6*t7);
	long double t45 = t43 + t44;
	long double t17 = t16*t6;
	long double t18 = t17 + t2;
	long double t23 = pow(r,3);
	long double t19 = 1/t18;
	long double t46 = 2 + r;
	long double t47 = r*t46;
	long double t48 = t47 + t6;
	long double t29 = 3*r;
	long double t20 = pow(r,4);
	long double t21 = -4*t20;
	long double t22 = -6*t2*t6;
	long double t24 = -2*t23*t6;
	long double t26 = -5*r*t25;
	long double t27 = 3 + r;
	long double t28 = t2*t27;
	long double t30 = 1 + t29;
	long double t31 = t30*t6;
	long double t32 = t28 + t31;
	long double t33 = 2*t16*t32*t6;
	long double t34 = -r;
	long double t35 = 1 + t34;
	long double t36 = 4*th;
	long double t37 = cos(t36);
	long double t38 = t25*t35*t37;
	long double t39 = t21 + t22 + t24 + t25 + t26 + t33 + t38;
	long double t49 = pow(t48,2);
	long double t52 = pow(t51,3);
	long double t53 = pow(t40,4);
	long double t54 = t2*t49*t52*t53*t6;
	long double t55 = pow(b1,2);
	long double t56 = pow(t18,2);
	long double t57 = t55*t56*t6*t7;
	long double t58 = r*t12*t6;
	long double t59 = t23 + t58;
	long double t60 = pow(t59,2);
	long double t61 = -4*t41*t6;
	long double t62 = t5 + t6 + t61;
	long double t63 = -(t60*t62);
	long double t64 = t57 + t63;
	long double t65 = t41*t45*t51*t64;
	long double t66 = t54 + t65;
	long double t67 = 1/t66;
	long double t69 = pow(t51,2);
	long double t93 = t16*t6*t7;
	long double t8 = b1*r*spin;
	long double t9 = 1/r;
	long double t10 = pow(spin,3);
	long double t13 = b1*t10*t12*t9;
	long double t14 = t13 + t8;
	long double t78 = 1 + r;
	long double t72 = 2*t2;
	long double t112 = t17 + t6 + t72;
	long double t146 = pow(r,9);
	long double t149 = pow(r,10);
	long double t152 = pow(r,6);
	long double t155 = pow(r,7);
	long double t144 = pow(r,8);
	long double t160 = pow(spin,6);
	long double t163 = pow(r,5);
	long double t168 = pow(spin,8);
	long double t174 = pow(spin,10);
	long double t182 = -15*r;
	long double t222 = 4*t55;
	long double t242 = 6*th;
	long double t243 = cos(t242);
	long double t177 = pow(spin,12);
	long double t255 = 8*th;
	long double t256 = cos(t255);
	long double t122 = pow(t40,3);
	long double t70 = r + t6;
	long double t71 = 2*t23*t70;
	long double t73 = t29 + t6 + t72;
	long double t74 = -2*r*t16*t6*t73;
	long double t75 = 2*t6;
	long double t76 = t2 + t29 + t75;
	long double t77 = -(r*t76);
	long double t79 = t16*t6*t78;
	long double t80 = t77 + t79;
	long double t81 = -2*t12*t6*t80;
	long double t82 = t71 + t74 + t81;
	long double t87 = 2*t43*t49;
	long double t88 = 6 + r;
	long double t89 = -(t23*t88);
	long double t90 = 4 + t29;
	long double t91 = -(r*t6*t90);
	long double t92 = -2*t25;
	long double t94 = t89 + t91 + t92 + t93;
	long double t95 = 2*t20;
	long double t96 = 2 + t29;
	long double t97 = r*t6*t96;
	long double t98 = t25 + t93 + t95 + t97;
	long double t99 = -(t94*t98);
	long double t100 = t87 + t99;
	long double t101 = sin(t15);
	long double t287 = pow(t18,-3);
	long double t103 = -6*t23*t6;
	long double t104 = -3*r*t25;
	long double t105 = 3*t2;
	long double t106 = t105 + t6;
	long double t107 = 2*t106*t16*t6*t78;
	long double t108 = t25*t37*t78;
	long double t109 = t103 + t104 + t107 + t108 + t21 + t22 + t25;
	long double t288 = -(t41*t49*t6);
	long double t289 = t45*t62;
	long double t290 = t288 + t289;
	long double t118 = t34 + t6;
	long double t124 = t23*t27;
	long double t125 = 2*r;
	long double t126 = 1 + t125;
	long double t127 = r*t126*t6;
	long double t128 = t124 + t127 + t25;
	long double t129 = r*t128;
	long double t130 = -t2;
	long double t131 = t130 + t6;
	long double t132 = -(t131*t16*t6);
	long double t133 = t129 + t132;
	long double t140 = pow(r,11);
	long double t142 = pow(r,12);
	long double t179 = -2*t2;
	long double t180 = -2 + r + t179;
	long double t181 = 32*t144*t180;
	long double t183 = 2 + t182;
	long double t184 = r*t183;
	long double t185 = -12*r;
	long double t186 = 5*t2;
	long double t187 = 4 + t185 + t186;
	long double t188 = t187*t55;
	long double t189 = t184 + t188;
	long double t190 = 8*t152*t189*t6;
	long double t191 = -64*t55;
	long double t192 = -12*t55;
	long double t193 = 7 + t192;
	long double t194 = 8*r*t193;
	long double t195 = 32*t55;
	long double t196 = -27 + t195;
	long double t197 = 3*t196*t2;
	long double t198 = 4 + t191 + t194 + t197;
	long double t199 = t198*t20*t25;
	long double t200 = 8 + t182;
	long double t201 = r*t200;
	long double t202 = -4*r;
	long double t203 = 13*t2;
	long double t204 = 4 + t202 + t203;
	long double t205 = 3*t204*t55;
	long double t206 = t201 + t205;
	long double t207 = 2*t160*t2*t206;
	long double t208 = -7*r;
	long double t209 = 6 + t208;
	long double t210 = -4*t209*t55;
	long double t211 = -5*r;
	long double t212 = t210 + t211;
	long double t213 = r*t168*t212;
	long double t214 = 6*t174*t55;
	long double t215 = t181 + t190 + t199 + t207 + t213 + t214;
	long double t217 = 10*r;
	long double t218 = -17*t2;
	long double t219 = -16 + t217 + t218;
	long double t220 = t152*t219;
	long double t221 = 8*t55;
	long double t223 = -15 + t222;
	long double t224 = t2*t223;
	long double t225 = 6*r*t55;
	long double t226 = r + t225;
	long double t227 = -2*t226;
	long double t228 = -2 + t221 + t224 + t227;
	long double t229 = 2*t20*t228*t6;
	long double t230 = 8*t2;
	long double t231 = 12*r;
	long double t232 = -9*t2;
	long double t233 = 4 + t231 + t232;
	long double t234 = t233*t55;
	long double t235 = t230 + t234;
	long double t236 = -2*t2*t235*t25;
	long double t237 = -1 + t222;
	long double t238 = 3*t160*t2*t237;
	long double t239 = 2*t168*t55;
	long double t240 = t220 + t229 + t236 + t238 + t239;
	long double t263 = 16*t155;
	long double t264 = -8*t20*t6;
	long double t265 = -18*t2*t25;
	long double t266 = -6*t23*t25;
	long double t267 = 6*r;
	long double t268 = 1 + t267;
	long double t269 = 8*t20*t268;
	long double t270 = 5*r;
	long double t271 = 3 + t270;
	long double t272 = 8*t2*t271*t6;
	long double t273 = 7*r;
	long double t274 = 1 + t273;
	long double t275 = t25*t274;
	long double t276 = t269 + t272 + t275;
	long double t277 = t16*t276*t6;
	long double t278 = -3*r;
	long double t279 = 7*t2;
	long double t280 = 4*t6;
	long double t281 = t278 + t279 + t280;
	long double t282 = 2*r*t25*t281*t37;
	long double t283 = -(t160*t243);
	long double t284 = r*t160*t243;
	long double t285 = t263 + t264 + t265 + t266 + t277 + t282 + t283 + t284;
	long double t361 = 1/t51;
	long double t111 = 1/t7;
	long double t113 = -3 + r;
	long double t139 = pow(t112,2);
	
	christ[0][0][0] = (t14*t19*t2*t39*t41*t45*t67*t7)/4.;
	christ[0][0][1] = -0.5*(t19*t2*t41*t67*t69*(-0.5*(t39*t45) + t41*t48*t6*t82));
	christ[0][0][2] = (t100*t101*t19*t2*t41*t6*t67*t69)/4.;
	christ[0][0][3] = (spin*t109*t14*t19*t2*t45*t53*t67*t7)/4.;
	christ[0][1][1] = (b1*spin*t111*t112*t41*t45*t56*(2*t113*t23 + t25 + r*(-1 + 4*r)*t6 + t118*t16*t6)*t67)/4.;
	christ[0][1][2] = -(t11*t122*t14*t2*t45*t56*t6*t67);
	christ[0][1][3] = -0.5*(spin*t112*t133*t53*t60*t67);
	christ[0][2][2] = -(t14*t23*t41*t45*t56*t67*t7);
	christ[0][2][3] = (8*t10*t11*t122*t139*t2*(-4*t2 + t20 + t25 + 2*t2*t6))/(64*t140 - 32*t142 - 24*t152*t160 + 144*t160*t163 - 3*t174*t2 + 16*t160*t20 - 15*t168*t20 + 30*t168*t23 - 14*t152*t160*t243 + 16*t160*t163*t243 - 6*t174*t2*t243 - 8*t160*t20*t243 - 20*t168*t20*t243 - 12*t144*t25 + 64*t152*t25 + 312*t155*t25 - t174*t2*t256 - t168*t20*t256 + 2*t168*t23*t256 + 4*t240*t25*t37 + 96*t152*t160*t55 - 96*t160*t163*t55 + 8*t177*t55 - 32*t168*t2*t55 + 48*t174*t2*t55 + 64*t160*t20*t55 + 88*t168*t20*t55 - 96*t168*t23*t55 - 16*r*t174*t243*t55 + 4*t177*t243*t55 + 16*t168*t2*t243*t55 + 8*t174*t2*t243*t55 + 4*t168*t20*t243*t55 - 16*t168*t23*t243*t55 + 80*t144*t25*t55 - 64*t152*t25*t55 - 64*t155*t25*t55 + 128*t144*t6 + 192*t146*t6 - 32*t149*t6 + 2*t16*t215*t6 - 64*t146*t55*t6 + 32*t149*t55*t6);
	christ[0][3][3] = -0.0625*(t14*t19*t2*t285*t45*t53*t67*t7);

	christ[1][0][0] = -0.25*(t2*t287*t290*t39*t41*t67*t69*t7);
	christ[1][0][1] = -0.5*(t14*t19*t2*t41*t67*t7*((t39*t45)/2. - t41*t48*t6*t82));
	christ[1][0][2] = -0.25*(t100*t101*t14*t19*t2*t41*t6*t67*t7);
	christ[1][0][3] = -0.25*(spin*t109*t287*t290*t53*t60*t67*t7);
	christ[1][1][1] = -0.5*(t2*t41*t51*t67*((2*t45*t55*t56*t6)/pow(r,3) - (t51*(2*r*t118 + 2*t12*t35*t6)*(t41*t49*t6 - t45*t62))/pow(t7,2))*t7);
	christ[1][1][2] = t11*t122*t2*t290*t6*t67*t69;
	christ[1][1][3] = (spin*t112*t133*t14*t2*t53*t67*t7)/2.;
	christ[1][2][2] = t23*t290*t41*t67*t69*t7;
	christ[1][2][3] = (32*b1*r*t11*t122*t25*t48*t51*pow(t7,2))/(-64*t140 + 32*t142 + 24*t152*t160 - 144*t160*t163 + 3*t174*t2 - 16*t160*t20 + 15*t168*t20 - 30*t168*t23 + 14*t152*t160*t243 - 16*t160*t163*t243 + 6*t174*t2*t243 + 8*t160*t20*t243 + 20*t168*t20*t243 + 12*t144*t25 - 64*t152*t25 - 312*t155*t25 + t174*t2*t256 + t168*t20*t256 - 2*t168*t23*t256 - 4*t240*t25*t37 - 96*t152*t160*t55 + 96*t160*t163*t55 - 8*t177*t55 + 32*t168*t2*t55 - 48*t174*t2*t55 - 64*t160*t20*t55 - 88*t168*t20*t55 + 96*t168*t23*t55 + 16*r*t174*t243*t55 - 4*t177*t243*t55 - 16*t168*t2*t243*t55 - 8*t174*t2*t243*t55 - 4*t168*t20*t243*t55 + 16*t168*t23*t243*t55 - 80*t144*t25*t55 + 64*t152*t25*t55 + 64*t155*t25*t55 - 128*t144*t6 - 192*t146*t6 + 32*t149*t6 - 2*t16*t215*t6 + 64*t146*t55*t6 - 32*t149*t55*t6);
	christ[1][3][3] = (t285*t287*t290*t53*t60*t67*t7)/16.;

	christ[2][0][0] = -0.5*(t101*t287*t361*t6*t94);
	christ[2][0][1] = 0;
	christ[2][0][2] = 0;
	christ[2][0][3] = (spin*t101*t287*t361*t43*t48)/2.;
	christ[2][1][1] = t11*t111*t361*t40*t6;
	christ[2][1][2] = r*t361;
	christ[2][1][3] = 0;
	christ[2][2][2] = -(t11*t361*t40*t6);
	christ[2][2][3] = 0;
	christ[2][3][3] = t287*t40*(-(t11*t18*t45) - 2*t101*t40*t45*t6 + t11*t18*t361*t41*t45*t6 + t11*t18*t41*t6*t7);
	
	christ[3][0][0] = (b1*r*t19*t39*t41*t48*t51*t6*t67*t7)/4.;
	christ[3][0][1] = -0.25*(spin*t19*t41*t67*(-(t39*t48*t60) - 2*t64*t82));
	christ[3][0][2] = (spin*t101*t48*t67*t7*(t16*(t160*(15 + t191)*t2 + 16*t20*t25*(3 - 2*t55) - 32*t168*t55 + 48*t152*t6) + t2*(32*t152 + 10*t160 + t160*t243 + 36*t2*t25 - 32*t160*t55 - 64*t2*t25*t55 + 48*t20*t6 - 32*t20*t55*t6 + 6*t25*t37*(t6 + t72))))/64.;
	christ[3][0][3] = (b1*r*t10*t109*t19*t48*t51*t53*t67*t7)/4.;
	christ[3][1][1] = (b1*t111*pow(t18,4)*t2*t41*t48*t51*t6*t67*(t118*t12*t6 + t2*(r*t113 + t75)))/pow(t23 + r*t16*t6,2);
	christ[3][1][2] = -(b1*r*t11*t122*t25*t48*t51*t56*t67);
	christ[3][1][3] = -0.5*(t19*t67*(-(t64*(2*r*t41*(t2*(t20 - t25 + t118*t41*t6) + t16*(t160 + 4*t2*t25 + 3*t20*t6 - t25*t41*(t278 + t6 + t72))) - 2*t12*t25*t53*(-(t16*t35*t6) + r*(t75 + 3*r*t78)))) + t2*t48*t53*t6*t69*t82));
	christ[3][2][2] = -(b1*t2*t41*t48*t51*t56*t6*t67*t7);
	christ[3][2][3] = t19*t67*((t11*t122*t139*t2*t43*t49*t6)/4. - t40*t64*(-(t11*t18*t45*t51) + t11*t18*t41*t45*t6 - 2*t101*t40*t45*t51*t6 + t11*t18*t41*t51*t6*t7));
	christ[3][3][3] = -0.0625*(b1*r*t19*t285*t48*t51*t53*t6*t67*t7);
};

void diffeqs(long double b, long double vars[], long double diffs[])
{
	long double t, r, th, phi;
	long double dt, dr, dth, dphi;
	long double christ[4][4][4];

	t = vars[0];
	r = vars[1];
	th = vars[2];
	phi = vars[3];

	dt = vars[4];
	dr = vars[5];
	dth = vars[6];
	dphi = vars[7];

	/* 1st order diff eqs */
	diffs[0] = dt; // dt
	diffs[1] = dr; // dr
	diffs[2] = dth; // dth
	diffs[3] = dphi; // dphi

	/* 2nd order diff eqs for r and theta; c.f. Eq (28) and Eq (29) */
	diffs[4] = -christ[0][0][0]*dt*dt - 2.0*christ[0][0][1]*dt*dr - 2.0*christ[0][0][2]*dt*dth - 2.0*christ[0][0][3]*dt*dphi 
			   -christ[0][1][1]*dr*dr - 2.0*christ[0][1][2]*dr*dth - 2.0*christ[0][1][3]*dr*dphi 
			   -christ[0][2][2]*dth*dth - 2.0*christ[0][2][3]*dth*dphi 
			   -christ[0][3][3]*dphi*dphi; // d2t
	diffs[5] = -christ[1][0][0]*dt*dt - 2.0*christ[1][0][1]*dt*dr - 2.0*christ[1][0][2]*dt*dth - 2.0*christ[1][0][3]*dt*dphi 
			   -christ[1][1][1]*dr*dr - 2.0*christ[1][1][2]*dr*dth - 2.0*christ[1][1][3]*dr*dphi 
			   -christ[1][2][2]*dth*dth - 2.0*christ[1][2][3]*dth*dphi 
			   -christ[1][3][3]*dphi*dphi; // d2r
	diffs[6] = -christ[2][0][0]*dt*dt - 2.0*christ[2][0][1]*dt*dr - 2.0*christ[2][0][2]*dt*dth - 2.0*christ[2][0][3]*dt*dphi 
			   -christ[2][1][1]*dr*dr - 2.0*christ[2][1][2]*dr*dth - 2.0*christ[2][1][3]*dr*dphi 
			   -christ[2][2][2]*dth*dth - 2.0*christ[2][2][3]*dth*dphi 
			   -christ[2][3][3]*dphi*dphi; // d2th
	diffs[7] = -christ[3][0][0]*dt*dt - 2.0*christ[3][0][1]*dt*dr - 2.0*christ[3][0][2]*dt*dth - 2.0*christ[3][0][3]*dt*dphi 
			   -christ[3][1][1]*dr*dr - 2.0*christ[3][1][2]*dr*dth - 2.0*christ[3][1][3]*dr*dphi 
			   -christ[3][2][2]*dth*dth - 2.0*christ[3][2][3]*dth*dphi 
			   -christ[3][3][3]*dphi*dphi; // d2phi
}
