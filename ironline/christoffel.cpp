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
double t13 = 2*t12;
double t18 = -2 + r;
double t19 = r*t18;
double t20 = t19 + t3;
double t8 = 2*t1;
double t23 = pow(r,6);
double t24 = 8*t23;
double t25 = 8*t12*t3;
double t27 = 8*t12*t26*t3;
double t28 = 8*r*t17*t26;
double t29 = 3*t1*t17;
double t30 = 12*t1*t17*t26;
double t31 = pow(spin,6);
double t32 = 4*t26*t31;
double t33 = t18*t26;
double t34 = r + t33;
double t35 = r*t3*t34;
double t36 = t17*t26;
double t37 = t13 + t35 + t36;
double t38 = 4*t3*t37*t5;
double t39 = 4*th;
double t40 = cos(t39);
double t41 = t1*t17*t40;
double t42 = t24 + t25 + t27 + t28 + t29 + t30 + t32 + t38 + t41;
double t43 = 1/t42;
double t2 = -2*t1;
double t7 = t2 + t3 + t6;
double t9 = t3 + t6 + t8;
double t11 = pow(t9,-2);
double t14 = 3*r;
double t15 = 2 + t14;
double t16 = r*t15*t3;
double t21 = t20*t3*t5;
double t22 = t13 + t16 + t17 + t21;
double t45 = 1/t20;
double t46 = t1 + t3;
double t47 = -t3;
double t48 = -(t3*t5);
double t49 = t47 + t48 + t8;
double t50 = 4*t1*t43*t45*t46*t49;
double t51 = pow(r,3);
double t52 = sin(t4);
double t54 = sin(th);
double t55 = pow(t54,2);
double t53 = -8*t3*t43*t51*t52;
double t68 = pow(spin,3);
double t69 = -2*defpar*r*t22*t43*t45*t52*t68;
double t56 = 8*defpar*r*t11*t22*t3*t43*t55*t7;
double t70 = -6*t12;
double t71 = -3*t1*t3;
double t72 = -t1;
double t73 = t3 + t72;
double t74 = t3*t5*t73;
double t75 = t17 + t70 + t71 + t74;
double t76 = 4*spin*t1*t43*t45*t55*t75;
double t78 = cos(th);
double t79 = pow(t54,3);
double t80 = 16*t43*t51*t68*t78*t79;
double t63 = -r;
double t98 = 1/t9;
double t100 = 8*defpar*r*spin*t43*t46*t7*t98;
double t81 = pow(r,5);
double t91 = 1 + t63;
double t101 = 16*defpar*t1*t20*t43*t52*t68*t98;
double t128 = -2*t1*t3*t43*t52*t9;
double t102 = 8*spin*t1*t20*t43*t55*t7*t98;
double t129 = -8*defpar*r*t3*t43*t55*t75*t98;
double t131 = -32*defpar*t1*t17*t20*t43*t78*t79*t98;
double t82 = 8*t81;
double t83 = -4*t1*t3;
double t84 = 8*t3*t51;
double t85 = 3*r*t17;
double t86 = 2*r;
double t87 = 1 + t86;
double t88 = r*t87;
double t89 = t3 + t88;
double t90 = 4*r*t3*t5*t89;
double t92 = -(t17*t40*t91);
double t93 = t17 + t82 + t83 + t84 + t85 + t90 + t92;
double t133 = pow(t9,-3);
double t136 = pow(t78,2);
double t137 = t136*t3;
double t138 = t1 + t137;
double t139 = 1/t138;
double t141 = r*t139;
double t135 = 8*r*spin*t133*t46*t52;
double t160 = 2*t17*t26;
double t161 = t1*t3*t5;
double t156 = -4 + t86;
double t157 = t156*t26;
double t158 = r + t157;
double t159 = r*t158*t3;
double t162 = t13 + t159 + t160 + t161;
double t163 = -4*spin*t162*t43*t45*t7*t98;
double t57 = pow(t20,-2);
double t58 = -3 + r;
double t59 = 2*t51*t58;
double t60 = 4*r;
double t61 = -1 + t60;
double t62 = r*t3*t61;
double t64 = t3 + t63;
double t65 = t3*t5*t64;
double t66 = t17 + t59 + t62 + t65;
double t164 = 2*t26;
double t165 = 1 + t164;
double t166 = t1*t165*t3;
double t167 = t13 + t160 + t161 + t166;
double t168 = 1/tan(th);
double t169 = -16*r*spin*t167*t168*t43*t98;
double t172 = -8*defpar*t1*t17*t43*t45*t52;
double t115 = pow(spin,8);
double t198 = t14 + t157;
double t170 = 32*defpar*t1*t11*t43*t55*t68*t7;
double t173 = 1/sin(th);
double t174 = pow(t173,2);
double t175 = 2*t18*t51;
double t176 = t159 + t160 + t161 + t175;
double t177 = -0.5*(t176*t55*t9*t93);
double t178 = pow(t54,4);
double t179 = -16*t138*t178*t3*t51*t7;
double t180 = t177 + t179;
double t181 = -2*t11*t174*t180*t43*t45;
double t183 = pow(r,8);
double t184 = 32*t183;
double t185 = 32*t3*t81;
double t186 = 48*t23*t3;
double t187 = 32*t23*t26*t3;
double t188 = 8*t17*t51;
double t189 = 64*t17*t26*t51;
double t190 = 36*t12*t17;
double t191 = 64*t12*t17*t26;
double t192 = 40*r*t26*t31;
double t193 = 10*t1*t31;
double t194 = 44*t1*t26*t31;
double t195 = 12*t115*t26;
double t196 = -2 + t14;
double t197 = 16*t196*t81;
double t199 = 16*t198*t3*t51;
double t200 = 15*r;
double t201 = 48*r;
double t202 = -32 + t201;
double t203 = t202*t26;
double t204 = t200 + t203;
double t205 = r*t17*t204;
double t206 = 16*t26*t31;
double t207 = t197 + t199 + t205 + t206;
double t208 = t207*t3*t5;
double t209 = -4*t51;
double t210 = 6*t12;
double t211 = r*t198*t3;
double t212 = t160 + t209 + t210 + t211;
double t213 = 2*t17*t212*t40;
double t214 = 6*th;
double t215 = cos(t214);
double t216 = t1*t215*t31;
double t217 = t184 + t185 + t186 + t187 + t188 + t189 + t190 + t191 + t192 + t193 + t194 + t195 + t208 + t213 + t216;
double t218 = (t168*t217*t43*t98)/2.;
christ[0][0][0] = -8*defpar*r*spin*t11*t22*t43*t7;
christ[0][0][1] = t50;
christ[0][0][2] = t53;
christ[0][0][3] = t56;
christ[0][1][0] = t50;
christ[0][1][1] = 2*defpar*spin*t22*t43*t57*t66;
christ[0][1][2] = t69;
christ[0][1][3] = t76;
christ[0][2][0] = t53;
christ[0][2][1] = t69;
christ[0][2][2] = -4*defpar*spin*t1*t22*t43;
christ[0][2][3] = t80;
christ[0][3][0] = t56;
christ[0][3][1] = t76;
christ[0][3][2] = t80;
christ[0][3][3] = -2*defpar*r*spin*t11*t22*t43*t55*t93;
christ[1][0][0] = -8*t1*t43*t49*(t47 + r*(2 + t63))*t98;
christ[1][0][1] = t100;
christ[1][0][2] = t101;
christ[1][0][3] = t102;
christ[1][1][0] = t100;
christ[1][1][1] = (t43*t45*(-8*pow(r,7) + t12*t17 - 4*t115*t26 + 16*t1*t17*t26 - 20*t12*t17*t26 + 4*t23*t3 - 8*t23*t26*t3 - 16*t1*t26*t31 + 3*t17*t51 + 16*t17*t26*t51 - 4*t3*t5*(t23 + 2*r*t17*t18*t26 + t26*t31 + t1*t3*(pow(t18,2)*t26 + t63)) + 16*t26*t3*t81 + t17*t40*t51*t91))/r;
christ[1][1][2] = t128;
christ[1][1][3] = t129;
christ[1][2][0] = t101;
christ[1][2][1] = t128;
christ[1][2][2] = -4*t20*t43*t51*t9;
christ[1][2][3] = t131;
christ[1][3][0] = t102;
christ[1][3][1] = t129;
christ[1][3][2] = t131;
christ[1][3][3] = -2*t1*t20*t43*t55*t93*t98;
christ[2][0][0] = -8*r*t133*t3*t52;
christ[2][0][1] = 0;
christ[2][0][2] = 0;
christ[2][0][3] = t135;
christ[2][1][0] = 0;
christ[2][1][1] = t139*t3*t45*t54*t78;
christ[2][1][2] = t141;
christ[2][1][3] = 0;
christ[2][2][0] = 0;
christ[2][2][1] = t141;
christ[2][2][2] = -(t139*t3*t54*t78);
christ[2][2][3] = 0;
christ[2][3][0] = t135;
christ[2][3][1] = 0;
christ[2][3][2] = 0;
christ[2][3][3] = -0.5*(t133*t52*(10*r*t17 + 11*t1*t17 + t24 + 16*t12*t3 + 3*t31 + t17*(-2*r + t1 + t3)*t40 + 16*t3*t51 + 4*t20*t3*t5*(t3 + t8)));
christ[3][0][0] = -32*defpar*t1*t11*t3*t43*t7;
christ[3][0][1] = t163;
christ[3][0][2] = t169;
christ[3][0][3] = t170;
christ[3][1][0] = t163;
christ[3][1][1] = 8*defpar*r*t3*t43*t57*t66;
christ[3][1][2] = t172;
christ[3][1][3] = t181;
christ[3][2][0] = t169;
christ[3][2][1] = t172;
christ[3][2][2] = -16*defpar*t3*t43*t51;
christ[3][2][3] = t218;
christ[3][3][0] = t170;
christ[3][3][1] = t181;
christ[3][3][2] = t218;
christ[3][3][3] = 8*defpar*t1*t11*t3*t43*t55*(-t17 - 3*r*t17 + 4*t1*t3 - 8*t3*t51 - 8*t81 - 4*r*t3*t5*t89 + t17*t40*t91);
}