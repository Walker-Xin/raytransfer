void redshift(long double r, long double th, long double ktkp, long double &gg)
{
long double t2 = pow(spin,2);
long double t1 = pow(r,2);
long double t13 = pow(spin,4);
long double t11 = pow(r,4);
long double t3 = 2*th;
long double t4 = cos(t3);
long double t35 = pow(spin,6);
long double t36 = 6*th;
long double t37 = cos(t36);
long double t5 = t2*t4;
long double t6 = t1 + t5;
long double t15 = pow(r,3);
long double t49 = 1 + r;
long double t32 = 4*th;
long double t33 = cos(t32);
long double t55 = sin(th);
long double t56 = pow(t55,2);
long double t8 = pow(r,7);
long double t9 = 16*t8;
long double t12 = -8*t11*t2;
long double t14 = -18*t1*t13;
long double t16 = -6*t13*t15;
long double t17 = 6*r;
long double t18 = 1 + t17;
long double t19 = 8*t11*t18;
long double t20 = 5*r;
long double t21 = 3 + t20;
long double t22 = 8*t1*t2*t21;
long double t23 = 7*r;
long double t24 = 1 + t23;
long double t25 = t13*t24;
long double t26 = t19 + t22 + t25;
long double t27 = t2*t26*t4;
long double t28 = -3*r;
long double t29 = 7*t1;
long double t30 = 4*t2;
long double t31 = t28 + t29 + t30;
long double t34 = 2*r*t13*t31*t33;
long double t38 = -(t35*t37);
long double t39 = r*t35*t37;
long double t40 = t12 + t14 + t16 + t27 + t34 + t38 + t39 + t9;
long double t45 = -4*t11;
long double t46 = -6*t1*t2;
long double t47 = -6*t15*t2;
long double t48 = -3*r*t13;
long double t50 = 3*t1;
long double t51 = t2 + t50;
long double t52 = 2*t2*t4*t49*t51;
long double t53 = t13*t33*t49;
long double t54 = t13 + t45 + t46 + t47 + t48 + t52 + t53;
long double t7 = pow(t6,3);
long double t41 = 1/t40;
long double t44 = pow(t6,-3);
long double t57 = 2*spin*t44*t54*t56;
long double t58 = pow(t6,-6);
long double t59 = 4*t11;
long double t60 = 6*t1*t2;
long double t61 = 2*t15*t2;
long double t62 = -t13;
long double t63 = 5*r*t13;
long double t64 = 3 + r;
long double t65 = t1*t64;
long double t66 = 3*r;
long double t67 = 1 + t66;
long double t68 = t2*t67;
long double t69 = t65 + t68;
long double t70 = -2*t2*t4*t69;
long double t71 = -r;
long double t72 = 1 + t71;
long double t73 = -(t13*t33*t72);
long double t74 = t59 + t60 + t61 + t62 + t63 + t70 + t73;
long double t75 = -(t40*t74);
long double t76 = pow(t54,2);
long double t77 = 4*t2*t56*t76;
long double t78 = t75 + t77;
long double t79 = t56*t58*t78;
long double t80 = sqrt(t79);
long double t81 = t57 + t80;
long double t90 = -2 + r;
long double t42 = 1/cos(th);
gg = sqrt(((-t2 + 4*t2*t56 + 4*spin*(r*(2 + r) + t2)*t41*t7*t81 - r*t90 - (4*pow(t6,6)*pow(2*spin*t44*t54*t55 + t42*t80,2)*(pow(t1 + t2,2) - t2*t56*(t2 + r*t90)))/pow(t40,2))*(t1 + t2*pow(cos(th),2)))/pow(t6,2))/(1 - 2*ktkp*t41*pow(t42,2)*t7*t81);
}

long double specific_energy(long double r)
{
long double t1 = pow(r,2);
long double t2 = -t1;
long double t3 = pow(spin,2);
long double t4 = t2 + t3;
long double t15 = pow(spin,4);
long double t28 = pow(r,6);
long double t23 = pow(r,8);
long double t33 = pow(r,4);
long double t9 = 3*r;
long double t37 = pow(spin,6);
long double t42 = pow(spin,8);
long double t45 = pow(spin,10);
long double t49 = pow(t4,-6);
long double t50 = pow(r,9);
long double t51 = -t50;
long double t52 = 2*t28;
long double t53 = -t23;
long double t54 = t52 + t53;
long double t55 = t3*t54;
long double t56 = 2 + t9;
long double t57 = 6*t15*t33*t56;
long double t58 = 16*r;
long double t59 = 10*t1;
long double t60 = 9 + t58 + t59;
long double t61 = 2*t1*t37*t60;
long double t62 = 4*r;
long double t63 = 5 + t62;
long double t64 = 3*r*t42*t63;
long double t65 = t45 + t51 + t55 + t57 + t61 + t64;
long double t66 = t1*t49*t65;
long double t67 = sqrt(t66);
long double t48 = 3 + r;
long double t44 = 1 + r;
long double t11 = 1 + t9;
long double t13 = 2*r;
long double t14 = 3 + t13;
long double se = (t1*(-((-2 + r)*r) + 3*t3 - (spin*pow(t1 - t3,3)*(r*(2 + r) + t3)*((8*r*spin*(pow(r,3) + t15 + 3*r*t3*t44))/pow(t4,3) + 8*t67))/(8.*(-pow(r,7) + t1*t14*t15 + t11*t3*t33))))/(pow(t4,2)*sqrt(-(((3 - r)*pow(r,11) - (29 + 40*r + 19*t1)*t15*t28 + (-9 - 8*r + 9*t1)*t23*t3 - (11 + 54*r + 35*t1)*t33*t37 - t1*(15 + 27*r + 16*t1)*t42 - 2*r*t44*t45 - 2*(9 + r)*pow(spin,3)*t23*t67 - 4*(-5 + r)*pow(spin,5)*t28*t67 + 4*(-3 + r)*pow(spin,7)*t33*t67 - 2*pow(spin,11)*t44*t67 + 2*pow(r,10)*spin*t48*t67 + 2*pow(spin,9)*t1*t48*t67)/(pow(-pow(r,5) + t14*t15 + t1*t11*t3,2)*t4))));
return se;
}