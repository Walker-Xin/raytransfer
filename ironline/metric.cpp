void metric(double r, double th, double g[4][4])
{
double t160 = pow(spin,2);
double t158 = pow(r,2);
double t167 = sin(th);
double t168 = pow(t167,2);
double t169 = -(t160*t168);
double t161 = cos(th);
double t162 = pow(t161,2);
double t180 = pow(r,3);
double t174 = t158 + t160;
double t202 = t160*t162;
double t204 = t158 + t202;
double t182 = defpar*t158;
double t197 = pow(r,4);
double t198 = 2*t197;
double t200 = defpar*t160;
double t201 = t182 + t198 + t200;
double t205 = defpar + t180;
double t206 = t174*t205;
double t207 = -(t160*t168*t180);
double t208 = t206 + t207;
double t209 = pow(t208,-2);
double t210 = -(spin*t168*t180*t201*t204*t209);
double t171 = pow(r,-3);
double t172 = defpar*t171;
double t173 = 1 + t172;
double t175 = t173*t174;
double t176 = t169 + t175;
double t178 = pow(t176,-2);
double t165 = -2 + r;
double t166 = r*t165;
g[0][0] = (-t158 - t160*t162)*(t160 + t166 + t169)*t178;
g[0][1] = 0;
g[0][2] = 0;
g[0][3] = t210;
g[1][0] = 0;
g[1][1] = t204/(-2*r + t158 + t160);
g[1][2] = 0;
g[1][3] = 0;
g[2][0] = 0;
g[2][1] = 0;
g[2][2] = t204;
g[2][3] = 0;
g[3][0] = t210;
g[3][1] = 0;
g[3][2] = 0;
g[3][3] = t168*(-(t160*(t160 + t166)*t168) + pow(t173,2)*pow(t174,2))*t178*t204;
}

void uppermetric(double r, double th, double gu[4][4])
{
double t159 = pow(spin,2);
double t158 = pow(r,2);
double t169 = -2 + r;
double t170 = r*t169;
double t171 = t159 + t170;
double t172 = 1/t171;
double t160 = cos(th);
double t161 = pow(t160,2);
double t162 = t159*t161;
double t163 = t158 + t162;
double t164 = 1/t163;
double t197 = pow(r,-3);
double t198 = defpar*t158;
double t200 = pow(r,4);
double t201 = 2*t200;
double t202 = defpar*t159;
double t204 = t198 + t201 + t202;
double t205 = -(spin*t164*t172*t197*t204);
double t176 = sin(th);
double t178 = pow(t176,2);
gu[0][0] = t164*(-((pow(defpar + pow(r,3),2)*pow(t158 + t159,2)*t172)/pow(r,6)) + t159*t178);
gu[0][1] = 0;
gu[0][2] = 0;
gu[0][3] = t205;
gu[1][0] = 0;
gu[1][1] = (-2*r + t158 + t159)*t164;
gu[1][2] = 0;
gu[1][3] = 0;
gu[2][0] = 0;
gu[2][1] = 0;
gu[2][2] = t164;
gu[2][3] = 0;
gu[3][0] = t205;
gu[3][1] = 0;
gu[3][2] = 0;
gu[3][3] = t164*t172*(t159 + t170 - t159*t178)*pow(1/sin(th),2);
}

void metric_rderivatives(double r, double th, double dg[4][4])
{
double t163 = pow(spin,2);
double t162 = pow(r,2);
double t166 = sin(th);
double t167 = pow(t166,2);
double t182 = -r;
double t176 = cos(th);
double t178 = pow(t176,2);
double t179 = t163*t178;
double t180 = t162 + t179;
double t164 = t162 + t163;
double t206 = pow(r,3);
double t168 = -(t163*t167);
double t207 = defpar + t206;
double t208 = t164*t207;
double t209 = -(t163*t167*t206);
double t210 = t208 + t209;
double t171 = defpar*t162;
double t219 = pow(r,4);
double t220 = 2*t219;
double t221 = defpar*t163;
double t222 = t171 + t220 + t221;
double t212 = -2 + r;
double t213 = r*t212;
double t205 = 1 + t182;
double t218 = pow(t210,-3);
double t223 = 2*defpar;
double t224 = 5*t206;
double t225 = 3*r*t163;
double t226 = -3*r*t163*t167;
double t227 = t223 + t224 + t225 + t226;
double t228 = 2*t162*t180*t222*t227;
double t229 = -2*t162*t210*t222;
double t230 = 2*defpar*r;
double t231 = 8*t206;
double t232 = t230 + t231;
double t234 = -(r*t180*t210*t232);
double t235 = -3*t180*t210*t222;
double t236 = t228 + t229 + t234 + t235;
double t237 = spin*t162*t167*t218*t236;
double t159 = pow(r,-3);
double t160 = defpar*t159;
double t161 = 1 + t160;
double t165 = t161*t164;
double t169 = t165 + t168;
double t170 = pow(t169,-3);
double t172 = pow(r,5);
double t173 = -2*t172;
double t174 = 3*defpar*t163;
double t175 = t171 + t173 + t174;
double t238 = t163 + t213;
double t256 = pow(t164,2);
dg[0][0] = (t170*(2*t175*t180*(-t163 + t163*t167 + r*(2 + t182)) + 2*r*t180*t205*t210 - 2*t162*t210*(t163 + t168 + t213)))/pow(r,4);
dg[0][1] = 0;
dg[0][2] = 0;
dg[0][3] = t237;
dg[1][0] = 0;
dg[1][1] = (2*(r*(t163 + t182) + t163*t178*t205))/pow(t238,2);
dg[1][2] = 0;
dg[1][3] = 0;
dg[2][0] = 0;
dg[2][1] = 0;
dg[2][2] = 2*r;
dg[2][3] = 0;
dg[3][0] = t237;
dg[3][1] = 0;
dg[3][2] = 0;
dg[3][3] = (2*t167*t170*(t180*(pow(r,7)*t163*t167*t205 - t164*t175*t207)*t210 + pow(r,11)*t169*(-(t163*t167*t238) + pow(t161,2)*t256) + t175*t180*(-(pow(r,6)*t163*t167*t238) + pow(t207,2)*t256)))/pow(r,10);
}

void metric_r2derivatives(double r, double th, double dg2[4][4])
{
double t159 = pow(defpar,2);
double t171 = pow(spin,2);
double t162 = pow(r,8);
double t166 = pow(r,11);
double t158 = pow(r,4);
double t182 = pow(spin,4);
double t173 = pow(r,6);
double t160 = pow(r,7);
double t176 = pow(r,9);
double t204 = pow(spin,6);
double t170 = pow(r,5);
double t180 = pow(r,3);
double t210 = pow(spin,8);
double t219 = 2*t173;
double t206 = pow(r,2);
double t239 = 4*defpar;
double t280 = 6*th;
double t281 = cos(t280);
double t164 = pow(r,10);
double t168 = pow(r,13);
double t420 = pow(defpar,3);
double t425 = pow(r,14);
double t433 = pow(r,12);
double t213 = -2*r;
double t214 = 1 + t213;
double t222 = -r;
double t246 = 2*th;
double t247 = cos(t246);
double t268 = 4*th;
double t276 = cos(t268);
double t414 = sin(th);
double t415 = pow(t414,2);
double t288 = defpar + t180;
double t412 = t171 + t206;
double t413 = t288*t412;
double t416 = -(t171*t180*t415);
double t417 = t413 + t416;
double t418 = pow(t417,-4);
double t421 = 96*t162*t420;
double t422 = 640*t159*t164;
double t423 = -576*t159*t166;
double t424 = -1024*defpar*t168;
double t426 = 192*defpar*t425;
double t427 = pow(r,16);
double t428 = 64*t427;
double t429 = 384*t171*t173*t420;
double t430 = 2272*t159*t162*t171;
double t431 = -2112*t159*t171*t176;
double t432 = -2688*defpar*t166*t171;
double t434 = 624*defpar*t171*t433;
double t435 = -64*t171*t425;
double t436 = 592*t158*t182*t420;
double t437 = 2688*t159*t173*t182;
double t438 = -2520*t159*t160*t182;
double t439 = -1152*defpar*t176*t182;
double t440 = 560*defpar*t164*t182;
double t441 = -72*t182*t433;
double t442 = 352*t204*t206*t420;
double t443 = 672*t158*t159*t204;
double t444 = -1104*t159*t170*t204;
double t445 = -144*defpar*t160*t204;
double t446 = 286*defpar*t162*t204;
double t447 = 48*t210*t420;
double t448 = -216*t159*t180*t210;
double t449 = 60*defpar*t173*t210;
double t450 = 6*t159*t214;
double t451 = -4*t173;
double t452 = 9*r;
double t453 = -16 + t452;
double t454 = defpar*t180*t453;
double t455 = t450 + t451 + t454;
double t456 = 16*t162*t455;
double t457 = 2 + t222;
double t458 = 36*defpar*t170*t457;
double t459 = 6*t162;
double t460 = 29*r;
double t461 = -12 + t460;
double t462 = 2*t159*t206*t461;
double t463 = t420 + t458 + t459 + t462;
double t464 = -16*t158*t171*t463;
double t465 = 32*t159;
double t466 = -12*r;
double t467 = 7 + t466;
double t468 = 96*defpar*t206*t467;
double t469 = 383*r;
double t470 = -192 + t469;
double t471 = t170*t470;
double t472 = t465 + t468 + t471;
double t473 = defpar*t182*t206*t472;
double t474 = 8*t159;
double t475 = -48*defpar*t180;
double t476 = 15*t173;
double t477 = t474 + t475 + t476;
double t478 = 6*defpar*t204*t477;
double t479 = t456 + t464 + t473 + t478;
double t480 = t171*t247*t479;
double t481 = -2*defpar*t180;
double t482 = 3*t170;
double t483 = t159 + t481 + t482;
double t484 = 4*t158*t483;
double t485 = 24*defpar;
double t486 = -49*r;
double t487 = 24 + t486;
double t488 = t206*t487;
double t489 = t485 + t488;
double t490 = defpar*t171*t206*t489;
double t491 = 2*defpar;
double t492 = -t180;
double t493 = t491 + t492;
double t494 = 18*defpar*t182*t493;
double t495 = t484 + t490 + t494;
double t496 = -2*t180*t182*t276*t495;
double t497 = defpar*t162*t204*t281;
double t498 = 6*defpar*t173*t210*t281;
double t499 = t421 + t422 + t423 + t424 + t426 + t428 + t429 + t430 + t431 + t432 + t434 + t435 + t436 + t437 + t438 + t439 + t440 + t441 + t442 + t443 + t444 + t445 + t446 + t447 + t448 + t449 + t480 + t496 + t497 + t498;
double t500 = -0.0625*(r*spin*t415*t418*t499);
double t515 = pow(defpar,4);
double t519 = pow(r,17);
double t528 = pow(r,15);
double t556 = pow(spin,10);
double t564 = 3*r;
double t507 = 6*r;
double t230 = -10*r;
dg2[0][0] = (t158*(320*t159*t160 - 240*t159*t162 - 512*defpar*t164 + 192*defpar*t166 + 32*t168 - 1344*defpar*t162*t171 - 32*t166*t171 + 1136*t159*t170*t171 - 864*t159*t171*t173 + 624*defpar*t171*t176 - 1078*t158*t159*t182 + 560*defpar*t160*t182 - 576*defpar*t173*t182 - 36*t176*t182 + 1344*t159*t180*t182 - 72*defpar*t158*t204 + 336*r*t159*t204 + 286*defpar*t170*t204 - 472*t159*t204*t206 - 90*t159*t210 + 60*defpar*t180*t210 - t171*(16*t170*(defpar*(8 - 9*r)*t180 - 3*t159*t214 + t219) + 24*t171*t180*((-8 + 15*r)*t159 + t219 + 24*defpar*t180*(1 + t222)) + defpar*r*t182*((96 - 383*r)*t180 - 48*defpar*(7 + t230)) + 30*defpar*t204*(-3*t180 + t239))*t247 - 2*t182*(t158*t159 - 8*defpar*t160 + 6*t176 + 3*defpar*(5*defpar - 6*t180)*t182 + defpar*t171*t206*(-49*t180 + 12*t206 + t239))*t276 + defpar*t170*t204*t281 + 6*defpar*t180*t210*t281)*t418)/8.;
dg2[0][1] = 0;
dg2[0][2] = 0;
dg2[0][3] = t500;
dg2[1][0] = 0;
dg2[1][1] = (2*(2*t180 + t182 - 3*t171*t206 - t171*(-4 + t171 - 3*t206 + t507)*pow(cos(th),2)))/pow((-2 + r)*r + t171,3);
dg2[1][2] = 0;
dg2[1][3] = 0;
dg2[2][0] = 0;
dg2[2][1] = 0;
dg2[2][2] = 2;
dg2[2][3] = 0;
dg2[3][0] = t500;
dg2[3][1] = 0;
dg2[3][2] = 0;
dg2[3][3] = (t415*t418*(128*pow(r,20) + 256*pow(r,18)*t171 + 1280*t159*t166*t171 - 2208*t159*t164*t182 + 1344*defpar*t168*t182 + 4448*t159*t176*t182 + 4992*t159*t160*t204 - 3288*t159*t162*t204 - 1152*defpar*t164*t204 + 736*defpar*t166*t204 - 48*t168*t204 - 96*defpar*t162*t210 + 672*t159*t170*t210 - 896*t159*t173*t210 + 160*defpar*t176*t210 + 12*t159*t162*t204*t281 + 16*defpar*t166*t204*t281 + 24*t168*t204*t281 + 48*defpar*t162*t210*t281 + 80*t159*t173*t210*t281 + 16*defpar*t176*t210*t281 + 512*t166*t420 + 2176*t171*t176*t420 + 3840*t160*t182*t420 + 3680*t170*t204*t420 + 1600*t180*t210*t420 + 768*t159*t425 - 2048*defpar*t171*t425 + 160*t204*t425 + 16*t204*t281*t425 + 288*t182*t427 + 960*t159*t171*t433 - 5120*defpar*t182*t433 + 35*t210*t433 + 8*t210*t281*t433 + 128*t162*t515 + 512*t171*t173*t515 + 768*t158*t182*t515 + 512*t204*t206*t515 + 128*t210*t515 + 512*defpar*t519 + 128*t171*t519 + 1280*defpar*t171*t528 - 64*t182*t528 - 168*t158*t159*t556 + 84*t158*t159*t281*t556 + 96*r*t420*t556 + 4*r*t182*t276*(-6*t159*t204*(-7*t180 + t239) + 2*t158*t171*(6*(1 + 2*r)*t162 + 36*defpar*(4 + r)*t170 + (-48 + 163*r)*t159*t206 + 4*t420) + t182*(-8*(21 - 40*r)*t158*t159 + 24*defpar*(1 + r)*t160 + 7*t166 - 16*t206*t420) + 8*t162*(-3*(1 - 3*r)*t159 + t173*(2 + t564) + 2*defpar*t180*(4 + t564))) - 4*t171*t180*t247*(21*r*t159*t210 + 2*t204*(6*defpar*(1 - 5*r)*t170 - 7*t176 - 38*t159*t180 + 128*t420) + t182*t206*(-252*defpar*t173 + (1152 - 1261*r)*t159*t206 + 6*t162*(1 + t230) + 416*t420) - 32*t158*t171*(3*t176 - 6*t420 + 2*defpar*t170*(19 + t507) + 17*t159*t206*(-2 + t564)) + 16*t173*((20 - 33*r)*t159*t206 + 2*t162*t214 + 2*t420 - 4*defpar*t170*(8 + t564))) + t210*t433*cos(8*th)))/64.;
}

void metric_alt(double spin, double spin2, double epsilon_r, double epsilon_t,
double z1, double z2, double mn[][4])
{
	double x, x2;
	double y, c2, s2, s4;
	double Delta;
	double Sigma, Sigma2;
	double h_t, h_r;
	double H,F;
	
	x = z1;
	y = z2;
	
	x2 = x*x;
	
	c2 = cos(y)*cos(y);
	s2 = 1 - c2;
	s4 = s2*s2;
	
	Delta = x2 - 2*x + spin2;
	
	Sigma  = x2 + spin2*c2;
	Sigma2 = Sigma*Sigma;
	
	h_t = epsilon_t*x/Sigma2;
	h_r = epsilon_r*x/Sigma2;
	
	H = (1+h_r)*(1+h_t);
	H = sqrt(H);
	F=1-2*x/Sigma;

	mn[0][0] = - (1 + h_t)*F;
	mn[0][3] = - spin*s2*(H-F*(1+h_t));
	mn[1][1] = Sigma*(1 + h_r)/(Delta + spin2*s2*h_r);
	mn[2][2] = Sigma;
	mn[3][0] = mn[0][3];
	mn[3][3] = s2*(Sigma+spin*spin*s2*(2*H-F*(1+h_t)));

	return ;
}