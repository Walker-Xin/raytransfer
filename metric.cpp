void metric(long double r, long double th, long double g[4][4])
{
long double t3 = pow(spin,2);
long double t1 = pow(r,2);
long double t11 = sin(th);
long double t12 = pow(t11,2);
long double t13 = -(t12*t3);
long double t4 = cos(th);
long double t5 = pow(t4,2);
long double t23 = pow(r,3);
long double t18 = t1 + t3;
long double t29 = t3*t5;
long double t30 = t1 + t29;
long double t24 = defpar*t1;
long double t25 = pow(r,4);
long double t26 = 2*t25;
long double t27 = defpar*t3;
long double t28 = t24 + t26 + t27;
long double t31 = defpar + t23;
long double t32 = t18*t31;
long double t33 = -(t12*t23*t3);
long double t34 = t32 + t33;
long double t35 = pow(t34,-2);
long double t36 = -(spin*t12*t23*t28*t30*t35);
long double t15 = pow(r,-3);
long double t16 = defpar*t15;
long double t17 = 1 + t16;
long double t19 = t17*t18;
long double t20 = t13 + t19;
long double t21 = pow(t20,-2);
long double t8 = -2 + r;
long double t9 = r*t8;
g[0][0] = t21*(-t1 - t3*t5)*(t13 + t3 + t9);
g[0][1] = 0;
g[0][2] = 0;
g[0][3] = t36;
g[1][0] = 0;
g[1][1] = t30/(-2*r + t1 + t3);
g[1][2] = 0;
g[1][3] = 0;
g[2][0] = 0;
g[2][1] = 0;
g[2][2] = t30;
g[2][3] = 0;
g[3][0] = t36;
g[3][1] = 0;
g[3][2] = 0;
g[3][3] = t12*t21*t30*(pow(t17,2)*pow(t18,2) - t12*t3*(t3 + t9));
}

void uppermetric(long double r, long double th, long double gu[4][4])
{
long double t2 = pow(spin,2);
long double t1 = pow(r,2);
long double t13 = -2 + r;
long double t14 = r*t13;
long double t15 = t14 + t2;
long double t16 = 1/t15;
long double t3 = cos(th);
long double t4 = pow(t3,2);
long double t5 = t2*t4;
long double t6 = t1 + t5;
long double t7 = 1/t6;
long double t25 = pow(r,-3);
long double t26 = defpar*t1;
long double t27 = pow(r,4);
long double t28 = 2*t27;
long double t29 = defpar*t2;
long double t30 = t26 + t28 + t29;
long double t31 = -(spin*t16*t25*t30*t7);
long double t20 = sin(th);
long double t21 = pow(t20,2);
gu[0][0] = (-((pow(defpar + pow(r,3),2)*t16*pow(t1 + t2,2))/pow(r,6)) + t2*t21)*t7;
gu[0][1] = 0;
gu[0][2] = 0;
gu[0][3] = t31;
gu[1][0] = 0;
gu[1][1] = (-2*r + t1 + t2)*t7;
gu[1][2] = 0;
gu[1][3] = 0;
gu[2][0] = 0;
gu[2][1] = 0;
gu[2][2] = t7;
gu[2][3] = 0;
gu[3][0] = t31;
gu[3][1] = 0;
gu[3][2] = 0;
gu[3][3] = t16*(t14 + t2 - t2*t21)*t7*pow(1/sin(th),2);
}

void metric_rderivatives(long double r, long double th, long double dg[4][4])
{
long double t6 = pow(spin,2);
long double t5 = pow(r,2);
long double t9 = sin(th);
long double t11 = pow(t9,2);
long double t24 = -r;
long double t20 = cos(th);
long double t21 = pow(t20,2);
long double t22 = t21*t6;
long double t23 = t22 + t5;
long double t7 = t5 + t6;
long double t32 = pow(r,3);
long double t12 = -(t11*t6);
long double t33 = defpar + t32;
long double t34 = t33*t7;
long double t35 = -(t11*t32*t6);
long double t36 = t34 + t35;
long double t15 = defpar*t5;
long double t45 = pow(r,4);
long double t46 = 2*t45;
long double t47 = defpar*t6;
long double t48 = t15 + t46 + t47;
long double t38 = -2 + r;
long double t39 = r*t38;
long double t31 = 1 + t24;
long double t44 = pow(t36,-3);
long double t49 = 2*defpar;
long double t50 = 5*t32;
long double t51 = 3*r*t6;
long double t52 = -3*r*t11*t6;
long double t53 = t49 + t50 + t51 + t52;
long double t54 = 2*t23*t48*t5*t53;
long double t55 = -2*t36*t48*t5;
long double t56 = 2*defpar*r;
long double t57 = 8*t32;
long double t58 = t56 + t57;
long double t59 = -(r*t23*t36*t58);
long double t60 = -3*t23*t36*t48;
long double t61 = t54 + t55 + t59 + t60;
long double t62 = spin*t11*t44*t5*t61;
long double t2 = pow(r,-3);
long double t3 = defpar*t2;
long double t4 = 1 + t3;
long double t8 = t4*t7;
long double t13 = t12 + t8;
long double t14 = pow(t13,-3);
long double t16 = pow(r,5);
long double t17 = -2*t16;
long double t18 = 3*defpar*t6;
long double t19 = t15 + t17 + t18;
long double t63 = t39 + t6;
long double t79 = pow(t7,2);
dg[0][0] = (t14*(2*r*t23*t31*t36 - 2*t36*t5*(t12 + t39 + t6) + 2*t19*t23*(r*(2 + t24) - t6 + t11*t6)))/pow(r,4);
dg[0][1] = 0;
dg[0][2] = 0;
dg[0][3] = t62;
dg[1][0] = 0;
dg[1][1] = (2*(t21*t31*t6 + r*(t24 + t6)))/pow(t63,2);
dg[1][2] = 0;
dg[1][3] = 0;
dg[2][0] = 0;
dg[2][1] = 0;
dg[2][2] = 2*r;
dg[2][3] = 0;
dg[3][0] = t62;
dg[3][1] = 0;
dg[3][2] = 0;
dg[3][3] = (2*t11*t14*(t23*t36*(pow(r,7)*t11*t31*t6 - t19*t33*t7) + t19*t23*(-(pow(r,6)*t11*t6*t63) + pow(t33,2)*t79) + pow(r,11)*t13*(-(t11*t6*t63) + pow(t4,2)*t79)))/pow(r,10);
}

void metric_r2derivatives(long double r, long double th, long double dg2[4][4])
{
long double t2 = pow(defpar,2);
long double t15 = pow(spin,2);
long double t5 = pow(r,8);
long double t9 = pow(r,11);
long double t1 = pow(r,4);
long double t24 = pow(spin,4);
long double t17 = pow(r,6);
long double t3 = pow(r,7);
long double t20 = pow(r,9);
long double t30 = pow(spin,6);
long double t14 = pow(r,5);
long double t23 = pow(r,3);
long double t36 = pow(spin,8);
long double t45 = 2*t17;
long double t32 = pow(r,2);
long double t64 = 4*defpar;
long double t87 = 6*th;
long double t88 = cos(t87);
long double t7 = pow(r,10);
long double t12 = pow(r,13);
long double t101 = pow(defpar,3);
long double t106 = pow(r,14);
long double t114 = pow(r,12);
long double t39 = -2*r;
long double t40 = 1 + t39;
long double t48 = -r;
long double t69 = 2*th;
long double t70 = cos(t69);
long double t84 = 4*th;
long double t85 = cos(t84);
long double t95 = sin(th);
long double t96 = pow(t95,2);
long double t92 = defpar + t23;
long double t93 = t15 + t32;
long double t94 = t92*t93;
long double t97 = -(t15*t23*t96);
long double t98 = t94 + t97;
long double t99 = pow(t98,-4);
long double t102 = 96*t101*t5;
long double t103 = 640*t2*t7;
long double t104 = -576*t2*t9;
long double t105 = -1024*defpar*t12;
long double t107 = 192*defpar*t106;
long double t108 = pow(r,16);
long double t109 = 64*t108;
long double t110 = 384*t101*t15*t17;
long double t111 = 2272*t15*t2*t5;
long double t112 = -2112*t15*t2*t20;
long double t113 = -2688*defpar*t15*t9;
long double t115 = 624*defpar*t114*t15;
long double t116 = -64*t106*t15;
long double t117 = 592*t1*t101*t24;
long double t118 = 2688*t17*t2*t24;
long double t119 = -2520*t2*t24*t3;
long double t120 = -1152*defpar*t20*t24;
long double t121 = 560*defpar*t24*t7;
long double t122 = -72*t114*t24;
long double t123 = 352*t101*t30*t32;
long double t124 = 672*t1*t2*t30;
long double t125 = -1104*t14*t2*t30;
long double t126 = -144*defpar*t3*t30;
long double t127 = 286*defpar*t30*t5;
long double t128 = 48*t101*t36;
long double t129 = -216*t2*t23*t36;
long double t130 = 60*defpar*t17*t36;
long double t131 = 6*t2*t40;
long double t132 = -4*t17;
long double t133 = 9*r;
long double t134 = -16 + t133;
long double t135 = defpar*t134*t23;
long double t136 = t131 + t132 + t135;
long double t137 = 16*t136*t5;
long double t138 = 2 + t48;
long double t139 = 36*defpar*t138*t14;
long double t140 = 6*t5;
long double t141 = 29*r;
long double t142 = -12 + t141;
long double t143 = 2*t142*t2*t32;
long double t144 = t101 + t139 + t140 + t143;
long double t145 = -16*t1*t144*t15;
long double t146 = 32*t2;
long double t147 = -12*r;
long double t148 = 7 + t147;
long double t149 = 96*defpar*t148*t32;
long double t150 = 383*r;
long double t151 = -192 + t150;
long double t152 = t14*t151;
long double t153 = t146 + t149 + t152;
long double t154 = defpar*t153*t24*t32;
long double t155 = 8*t2;
long double t156 = -48*defpar*t23;
long double t157 = 15*t17;
long double t158 = t155 + t156 + t157;
long double t159 = 6*defpar*t158*t30;
long double t160 = t137 + t145 + t154 + t159;
long double t161 = t15*t160*t70;
long double t162 = -2*defpar*t23;
long double t163 = 3*t14;
long double t164 = t162 + t163 + t2;
long double t165 = 4*t1*t164;
long double t166 = 24*defpar;
long double t167 = -49*r;
long double t168 = 24 + t167;
long double t169 = t168*t32;
long double t170 = t166 + t169;
long double t171 = defpar*t15*t170*t32;
long double t172 = 2*defpar;
long double t173 = -t23;
long double t174 = t172 + t173;
long double t175 = 18*defpar*t174*t24;
long double t176 = t165 + t171 + t175;
long double t177 = -2*t176*t23*t24*t85;
long double t178 = defpar*t30*t5*t88;
long double t179 = 6*defpar*t17*t36*t88;
long double t180 = t102 + t103 + t104 + t105 + t107 + t109 + t110 + t111 + t112 + t113 + t115 + t116 + t117 + t118 + t119 + t120 + t121 + t122 + t123 + t124 + t125 + t126 + t127 + t128 + t129 + t130 + t161 + t177 + t178 + t179;
long double t181 = -0.0625*(r*spin*t180*t96*t99);
long double t196 = pow(defpar,4);
long double t200 = pow(r,17);
long double t209 = pow(r,15);
long double t237 = pow(spin,10);
long double t245 = 3*r;
long double t188 = 6*r;
long double t56 = -10*r;
dg2[0][0] = (t1*(32*t12 + 1136*t14*t15*t2 - 864*t15*t17*t2 + 624*defpar*t15*t20 - 576*defpar*t17*t24 - 1078*t1*t2*t24 - 36*t20*t24 + 1344*t2*t23*t24 + 320*t2*t3 + 560*defpar*t24*t3 - 72*defpar*t1*t30 + 286*defpar*t14*t30 + 336*r*t2*t30 - 472*t2*t30*t32 - 90*t2*t36 + 60*defpar*t23*t36 - 1344*defpar*t15*t5 - 240*t2*t5 - 512*defpar*t7 - t15*(16*t14*(defpar*(8 - 9*r)*t23 - 3*t2*t40 + t45) + 24*t15*t23*((-8 + 15*r)*t2 + t45 + 24*defpar*t23*(1 + t48)) + defpar*r*t24*((96 - 383*r)*t23 - 48*defpar*(7 + t56)) + 30*defpar*t30*(-3*t23 + t64))*t70 - 2*t24*(t1*t2 + 6*t20 + 3*defpar*(5*defpar - 6*t23)*t24 - 8*defpar*t3 + defpar*t15*t32*(-49*t23 + 12*t32 + t64))*t85 + defpar*t14*t30*t88 + 6*defpar*t23*t36*t88 + 192*defpar*t9 - 32*t15*t9)*t99)/8.;
dg2[0][1] = 0;
dg2[0][2] = 0;
dg2[0][3] = t181;
dg2[1][0] = 0;
dg2[1][1] = (2*(2*t23 + t24 - 3*t15*t32 - t15*(-4 + t15 + t188 - 3*t32)*pow(cos(th),2)))/pow((-2 + r)*r + t15,3);
dg2[1][2] = 0;
dg2[1][3] = 0;
dg2[2][0] = 0;
dg2[2][1] = 0;
dg2[2][2] = 2;
dg2[2][3] = 0;
dg2[3][0] = t181;
dg2[3][1] = 0;
dg2[3][2] = 0;
dg2[3][3] = (t96*t99*(128*pow(r,20) + 256*pow(r,18)*t15 - 2048*defpar*t106*t15 + 512*t15*t17*t196 + 768*t106*t2 + 960*t114*t15*t2 + 2176*t101*t15*t20 + 512*defpar*t200 + 128*t15*t200 + 1280*defpar*t15*t209 + 96*r*t101*t237 - 168*t1*t2*t237 + 288*t108*t24 - 5120*defpar*t114*t24 + 1344*defpar*t12*t24 + 768*t1*t196*t24 + 4448*t2*t20*t24 - 64*t209*t24 + 3840*t101*t24*t3 + 160*t106*t30 - 48*t12*t30 + 3680*t101*t14*t30 + 4992*t2*t3*t30 + 512*t196*t30*t32 + 35*t114*t36 + 128*t196*t36 + 672*t14*t2*t36 - 896*t17*t2*t36 + 160*defpar*t20*t36 + 1600*t101*t23*t36 + 128*t196*t5 - 3288*t2*t30*t5 - 96*defpar*t36*t5 - 2208*t2*t24*t7 - 1152*defpar*t30*t7 - 4*t15*t23*(2*(128*t101 + 6*defpar*(1 - 5*r)*t14 - 7*t20 - 38*t2*t23)*t30 - 32*t1*t15*(-6*t101 + 2*defpar*t14*(19 + t188) + 3*t20 + 17*t2*(-2 + t245)*t32) + 21*r*t2*t36 + 16*t17*(2*t101 - 4*defpar*t14*(8 + t245) + (20 - 33*r)*t2*t32 + 2*t40*t5) + t24*t32*(416*t101 - 252*defpar*t17 + (1152 - 1261*r)*t2*t32 + 6*t5*(1 + t56)))*t70 + 84*t1*t2*t237*t88 + 16*t106*t30*t88 + 24*t12*t30*t88 + 8*t114*t36*t88 + 80*t17*t2*t36*t88 + 16*defpar*t20*t36*t88 + 12*t2*t30*t5*t88 + 48*defpar*t36*t5*t88 + 512*t101*t9 + 1280*t15*t2*t9 + 736*defpar*t30*t9 + 16*defpar*t30*t88*t9 + 4*r*t24*t85*(8*(-3*(1 - 3*r)*t2 + t17*(2 + t245) + 2*defpar*t23*(4 + t245))*t5 + 2*t1*t15*(4*t101 + 36*defpar*(4 + r)*t14 + (-48 + 163*r)*t2*t32 + 6*(1 + 2*r)*t5) - 6*t2*t30*(-7*t23 + t64) + t24*(-8*(21 - 40*r)*t1*t2 + 24*defpar*(1 + r)*t3 - 16*t101*t32 + 7*t9)) + t114*t36*cos(8*th)))/64.;
}