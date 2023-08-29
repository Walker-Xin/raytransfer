void christoffel(double spin, double defpar, double r, double th, double christ[4][4][4])
{
double t161 = pow(spin,2);
double t170 = pow(defpar,2);
double t173 = pow(r,8);
double t178 = pow(r,11);
double t208 = pow(spin,4);
double t200 = pow(r,6);
double t171 = pow(r,7);
double t204 = pow(r,9);
double t164 = pow(r,2);
double t215 = pow(spin,6);
double t210 = pow(r,4);
double t197 = pow(r,5);
double t207 = pow(r,3);
double t221 = pow(spin,8);
double t159 = -2 + r;
double t246 = -2*t171;
double t252 = defpar + t207;
double t165 = cos(th);
double t166 = pow(t165,2);
double t167 = t161*t166;
double t168 = t164 + t167;
double t169 = 1/t168;
double t258 = t170*t210;
double t255 = 2*th;
double t256 = cos(t255);
double t418 = sin(th);
double t416 = t161 + t164;
double t417 = t252*t416;
double t419 = pow(t418,2);
double t420 = -(t161*t207*t419);
double t421 = t417 + t420;
double t422 = pow(t421,-2);
double t158 = 1/r;
double t160 = r*t159;
double t162 = t160 + t161;
double t163 = 1/t162;
double t172 = -80*t170*t171;
double t174 = 48*t170*t173;
double t175 = pow(r,10);
double t176 = -64*defpar*t175;
double t179 = 48*defpar*t178;
double t180 = pow(r,13);
double t182 = 16*t180;
double t198 = -248*t161*t170*t197;
double t201 = 172*t161*t170*t200;
double t202 = -208*defpar*t161*t173;
double t205 = 160*defpar*t161*t204;
double t206 = 16*t161*t178;
double t209 = -224*t170*t207*t208;
double t211 = 215*t170*t208*t210;
double t212 = -188*defpar*t200*t208;
double t213 = 182*defpar*t171*t208;
double t214 = -6*t204*t208;
double t216 = -56*r*t170*t215;
double t217 = 112*t164*t170*t215;
double t218 = -56*defpar*t210*t215;
double t219 = 88*defpar*t197*t215;
double t220 = -6*t171*t215;
double t222 = 21*t170*t221;
double t223 = 18*defpar*t207*t221;
double t224 = 2*r;
double t225 = -3 + t224;
double t226 = 4*t207*t225;
double t227 = 5*r;
double t228 = -6 + t227;
double t229 = defpar*t228;
double t230 = t226 + t229;
double t231 = defpar*t197*t230;
double t232 = -9*r;
double t234 = 10 + t232;
double t235 = t170*t207*t234;
double t236 = -13*r;
double t237 = 16 + t236;
double t238 = defpar*t200*t237;
double t239 = t204 + t235 + t238;
double t240 = -2*t161*t239;
double t241 = -3*r;
double t242 = 2 + t241;
double t244 = -8*defpar*t210*t242;
double t247 = 19*r;
double t248 = -14 + t247;
double t249 = r*t170*t248;
double t250 = t244 + t246 + t249;
double t251 = t208*t250;
double t253 = 6*defpar*t215*t252;
double t254 = t231 + t240 + t251 + t253;
double t257 = 4*t161*t254*t256;
double t259 = 2*defpar*t159*t200;
double t260 = -2*t204;
double t268 = 4*t164*t170;
double t276 = -1 + r;
double t279 = 8*defpar*t210*t276;
double t280 = t246 + t268 + t279;
double t281 = t161*t280;
double t285 = 2*t207;
double t286 = defpar + t285;
double t287 = 3*defpar*t208*t286;
double t288 = t258 + t259 + t260 + t281 + t287;
double t412 = 4*th;
double t413 = cos(t412);
double t414 = t208*t288*t413;
double t415 = t172 + t174 + t176 + t179 + t182 + t198 + t201 + t202 + t205 + t206 + t209 + t211 + t212 + t213 + t214 + t216 + t217 + t218 + t219 + t220 + t222 + t223 + t257 + t414;
double t423 = (t158*t163*t169*t415*t422)/16.;
double t438 = pow(defpar,3);
double t439 = -2*r;
double t424 = 2*defpar*t171;
double t425 = 2*t204;
double t426 = 2*t164*t170;
double t427 = 3*defpar*t197;
double t428 = t171 + t426 + t427;
double t429 = t161*t428;
double t430 = defpar*t208*t252;
double t431 = defpar*t164;
double t432 = defpar*t161;
double t433 = t210 + t431 + t432;
double t434 = t161*t207*t256*t433;
double t435 = t258 + t424 + t425 + t429 + t430 + t434;
double t436 = -(t161*t165*t169*t418*t422*t435);
double t437 = pow(r,-4);
double t440 = 3 + t439;
double t441 = 2*t164*t170*t440;
double t442 = -5*defpar*t200;
double t443 = -6*t173;
double t444 = t438 + t441 + t442 + t443;
double t445 = t210*t444;
double t446 = 2*t438;
double t447 = 7 + t439;
double t448 = t164*t170*t447;
double t449 = -r;
double t450 = 2 + t449;
double t451 = 4*defpar*t197*t450;
double t452 = t173 + t446 + t448 + t451;
double t453 = 2*t161*t164*t452;
double t454 = -(defpar*t200);
double t455 = t438 + t454;
double t456 = 3*t208*t455;
double t457 = t445 + t453 + t456;
double t458 = t416*t457;
double t459 = 2*t200;
double t460 = 3*r;
double t461 = -4 + t460;
double t462 = defpar*t207*t461;
double t463 = 4 + t241;
double t464 = defpar*t463;
double t465 = t207 + t464;
double t466 = -2*r*t161*t465;
double t467 = 3*defpar*t208;
double t468 = t459 + t462 + t466 + t467;
double t469 = t161*t200*t419*t468;
double t470 = t458 + t469;
double t471 = (spin*t163*t419*t422*t437*t470)/2.;
double t472 = pow(spin,3);
double t473 = 2*t210;
double t474 = t431 + t432 + t473;
double t475 = pow(t418,3);
double t476 = t165*t207*t422*t472*t474*t475;
double t500 = pow(t421,-3);
double t522 = -(t161*t165*t169*t418);
double t502 = 2*defpar;
double t503 = 5*t207;
double t504 = 3*r*t161;
double t505 = -3*r*t161*t419;
double t506 = t502 + t503 + t504 + t505;
double t507 = -2*t164*t168*t474*t506;
double t508 = 2*t164*t421*t474;
double t509 = 2*defpar*r;
double t510 = 8*t207;
double t511 = t509 + t510;
double t512 = r*t168*t421*t511;
double t513 = 3*t168*t421*t474;
double t514 = t507 + t508 + t512 + t513;
double t515 = (spin*t162*t164*t169*t419*t500*t514)/2.;
double t518 = 1 + t449;
double t564 = sin(t255);
double t588 = r*t169;
double t580 = r*t161;
double t581 = t285 + t502 + t580;
double t582 = t164*t581;
double t583 = t207 + t502;
double t584 = t161*t256*t583;
double t585 = t582 + t584;
double t586 = (spin*t169*t207*t416*t474*t500*t564*t585)/4.;
double t589 = pow(r,-3);
double t590 = defpar*t589;
double t591 = 1 + t590;
double t562 = pow(t418,4);
double t592 = t416*t591;
double t593 = -(t161*t419);
double t594 = t592 + t593;
double t596 = pow(t591,2);
double t597 = pow(t416,2);
double t601 = t596*t597;
double t602 = -(t161*t162*t419);
double t603 = t601 + t602;
double t496 = 3*defpar*t161;
double t608 = -8*defpar*t207;
double t609 = 6*defpar*t210;
double t610 = 4*t200;
double t611 = -2*t207;
double t612 = 11*r;
double t613 = -16 + t612;
double t614 = defpar*t613;
double t615 = t611 + t614;
double t616 = r*t161*t615;
double t617 = -2*t210;
double t618 = t431 + t496 + t617;
double t619 = t161*t256*t618;
double t620 = t467 + t608 + t609 + t610 + t616 + t619;
double t621 = (spin*t163*t164*t422*t620)/4.;
double t552 = 3*t170*t208;
double t682 = 6*th;
double t683 = cos(t682);
double t622 = t197 + t431;
double t623 = 2*t622;
double t624 = t161*t583;
double t625 = t161*t207*t256;
double t626 = t623 + t624 + t625;
double t627 = pow(t626,-2);
double t628 = 1/tan(th);
double t629 = -4*spin*t207*t474*t627*t628;
double t665 = 32*defpar*t207;
double t630 = -64*t170*t171;
double t631 = 32*t170*t173;
double t632 = -128*defpar*t175;
double t633 = 64*defpar*t178;
double t634 = -64*t180;
double t635 = pow(r,14);
double t636 = 32*t635;
double t637 = -128*t161*t170*t197;
double t638 = 104*t161*t170*t200;
double t639 = -160*defpar*t161*t173;
double t640 = 128*defpar*t161*t204;
double t641 = -80*t161*t178;
double t642 = pow(r,12);
double t643 = 48*t161*t642;
double t644 = -64*t170*t207*t208;
double t645 = 130*t170*t208*t210;
double t646 = 8*defpar*t200*t208;
double t647 = 88*defpar*t171*t208;
double t648 = -24*t204*t208;
double t649 = 36*t175*t208;
double t650 = 64*t164*t170*t215;
double t651 = 16*defpar*t210*t215;
double t652 = 24*defpar*t197*t215;
double t653 = 2*t171*t215;
double t654 = 10*t173*t215;
double t655 = 6*t170*t221;
double t656 = 4*defpar*t164*t440;
double t657 = 6*t197*t518;
double t658 = t170 + t656 + t657;
double t659 = -8*t210*t658;
double t660 = 2*defpar*t210*t463;
double t661 = t171*t242;
double t662 = t426 + t660 + t661;
double t663 = -16*t161*t662;
double t664 = -24*t170;
double t666 = 15*r;
double t667 = 1 + t666;
double t668 = t197*t667;
double t669 = t664 + t665 + t668;
double t670 = t208*t669;
double t671 = t659 + t663 + t670;
double t672 = t161*t164*t256*t671;
double t673 = 4*defpar*t200*t518;
double t674 = 2*t204*t242;
double t675 = 4*defpar*t210*t450;
double t676 = 1 + t241;
double t677 = t171*t676;
double t678 = t268 + t675 + t677;
double t679 = t161*t678;
double t680 = t258 + t552 + t673 + t674 + t679;
double t681 = -2*t208*t413*t680;
double t684 = -(t171*t215*t683);
double t685 = t173*t215*t683;
double t686 = t630 + t631 + t632 + t633 + t634 + t636 + t637 + t638 + t639 + t640 + t641 + t643 + t644 + t645 + t646 + t647 + t648 + t649 + t650 + t651 + t652 + t653 + t654 + t655 + t672 + t681 + t684 + t685;
double t687 = (t158*t163*t169*t422*t686)/32.;
double t688 = 2*t170;
double t689 = 4*defpar*t207;
double t690 = -2*t197;
double t691 = 3*t200;
double t692 = t688 + t689 + t690 + t691;
double t693 = 16*t210*t692;
double t694 = 6*defpar*t197;
double t695 = 3*t173;
double t696 = t268 + t694 + t695;
double t697 = 16*t161*t696;
double t698 = 32*t170;
double t699 = 15*t200;
double t700 = t665 + t698 + t699;
double t701 = t208*t700;
double t702 = t693 + t697 + t701;
double t703 = t161*t256*t702;
double t704 = 32*t170*t210;
double t705 = 64*defpar*t171;
double t706 = 32*t175;
double t707 = 64*t161*t164*t170;
double t708 = 128*defpar*t161*t197;
double t709 = 32*t161*t171;
double t710 = 48*t161*t173;
double t711 = 32*t170*t208;
double t712 = 88*defpar*t207*t208;
double t713 = 8*t197*t208;
double t714 = 36*t200*t208;
double t715 = 24*defpar*r*t215;
double t716 = 10*t210*t215;
double t717 = 4*defpar*t164;
double t718 = -4*t210;
double t719 = 6*t197;
double t720 = 4*defpar;
double t721 = 3*t207;
double t722 = t720 + t721;
double t723 = t161*t722;
double t724 = t717 + t718 + t719 + t723;
double t725 = 2*r*t208*t413*t724;
double t726 = t210*t215*t683;
double t727 = t704 + t705 + t706 + t707 + t708 + t709 + t710 + t711 + t712 + t713 + t714 + t715 + t716 + t725 + t726;
double t728 = t164*t727;
double t729 = t703 + t728;
double t730 = (t169*t627*t628*t729)/8.;
christ[0][0][0] = 0;
christ[0][0][1] = t423;
christ[0][0][2] = t436;
christ[0][0][3] = 0;
christ[0][1][0] = t423;
christ[0][1][1] = 0;
christ[0][1][2] = 0;
christ[0][1][3] = t471;
christ[0][2][0] = t436;
christ[0][2][1] = 0;
christ[0][2][2] = 0;
christ[0][2][3] = t476;
christ[0][3][0] = 0;
christ[0][3][1] = t471;
christ[0][3][2] = t476;
christ[0][3][3] = 0;
christ[1][0][0] = (t162*t169*t197*(8*t173 - 40*defpar*t197 + 24*defpar*t200 - 84*defpar*t161*t207 - 28*defpar*r*t208 + 35*defpar*t164*t208 + 56*defpar*t161*t210 - 3*t208*t210 + 9*defpar*t215 + 4*t161*t256*(defpar*(-3 + 4*r)*t207 - r*t161*(t207 + defpar*(7 + t232)) + t467) + t208*t413*(-t210 + t431 + t496))*t500)/8.;
christ[1][0][1] = 0;
christ[1][0][2] = 0;
christ[1][0][3] = t515;
christ[1][1][0] = 0;
christ[1][1][1] = t163*t169*(r*(t161 + t449) + t161*t166*t518);
christ[1][1][2] = t522;
christ[1][1][3] = 0;
christ[1][2][0] = 0;
christ[1][2][1] = t522;
christ[1][2][2] = -(r*t162*t169);
christ[1][2][3] = 0;
christ[1][3][0] = t515;
christ[1][3][1] = 0;
christ[1][3][2] = 0;
christ[1][3][3] = -(r*t162*t169*t419*t500*(pow(t252,3)*pow(t416,3) + t419*(-(t161*t207*(defpar*(-5 + 6*r)*t200 + t161*(-2*t164*t170 + t171*(-3 + t227) + 9*defpar*t210*t276) + t208*(-2*t170 + 3*defpar*t207 + t459) + t204*(1 + t460))) + r*t166*t208*(-((3 + r)*t204) + t258 + 3*defpar*t200*t518 + t161*t210*(defpar*(7 + t241) + t207*t518) + t552)) + t173*t208*(r*(t161 + r*t225) - t161*t166*t518)*t562 + t170*t207*t215*pow(t564,2)));
christ[2][0][0] = -0.5*(t161*t169*t200*(defpar*t208 + 2*t207*(t207 + defpar*t276) + t161*t256*t433 + r*t161*(t207 + defpar*(-2 + t460)))*t500*t564);
christ[2][0][1] = 0;
christ[2][0][2] = 0;
christ[2][0][3] = t586;
christ[2][1][0] = 0;
christ[2][1][1] = t161*t163*t165*t169*t418;
christ[2][1][2] = t588;
christ[2][1][3] = 0;
christ[2][2][0] = 0;
christ[2][2][1] = t588;
christ[2][2][2] = t522;
christ[2][2][3] = 0;
christ[2][3][0] = t586;
christ[2][3][1] = 0;
christ[2][3][2] = 0;
christ[2][3][3] = (t165*t418*(2*t162*t208*t562 + t161*t162*t419*t594 - 2*t161*t419*t596*t597 - t594*t603 + t161*t169*t419*t594*t603))/pow(t594,3);
christ[3][0][0] = 0;
christ[3][0][1] = t621;
christ[3][0][2] = t629;
christ[3][0][3] = 0;
christ[3][1][0] = t621;
christ[3][1][1] = 0;
christ[3][1][2] = 0;
christ[3][1][3] = t687;
christ[3][2][0] = t629;
christ[3][2][1] = 0;
christ[3][2][2] = 0;
christ[3][2][3] = t730;
christ[3][3][0] = 0;
christ[3][3][1] = t687;
christ[3][3][2] = t730;
christ[3][3][3] = 0;
}

void christoffel_alt(double spin, double spin2, double epsilon_r, double epsilon_t,
double w1, double w2, double CS[][4][4])
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
	
	// metric_alt(spin, spin2, epsilon_r, epsilon_t, w1,w2,g);
	metric(w1,w2,g);
	
	gg = g[0][0]*g[3][3] - g[0][3]*g[0][3];
	
	invg[0][0] = g[3][3]/gg; 
	invg[0][3] = - g[0][3]/gg;
	invg[1][1] = 1/g[1][1];
	invg[2][2] = 1/g[2][2];
	invg[3][0] = invg[0][3]; 
	invg[3][3] = g[0][0]/gg;
	
	
	/* ----- compute metric derivatives ----- */
	
	
	// metric_alt(spin, spin2, epsilon_r, epsilon_t, w1+dw1,w2,gr);
	// metric_alt(spin, spin2, epsilon_r, epsilon_t, w1-dw1,w2,gl);
	metric(w1+dw1,w2,gr);
	metric(w1-dw1,w2,gl);
	
	Dg[0][0][1] = 0.5*(gr[0][0] - gl[0][0])/dw1;
	Dg[0][3][1] = 0.5*(gr[0][3] - gl[0][3])/dw1;
	Dg[1][1][1] = 0.5*(gr[1][1] - gl[1][1])/dw1;
	Dg[2][2][1] = 0.5*(gr[2][2] - gl[2][2])/dw1;
	Dg[3][0][1] = Dg[0][3][1];
	Dg[3][3][1] = 0.5*(gr[3][3] - gl[3][3])/dw1;
	
	// metric_alt(spin, spin2, epsilon_r, epsilon_t, w1,w2+dw2,gr);
	// metric_alt(spin, spin2, epsilon_r, epsilon_t, w1,w2-dw2,gl);
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
	
	CS[0][0][0] = 0;
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
	CS[3][3][3] = 0;

	return ;
}

void diffeqs_alt(double christ[4][4][4], double vars[], double diffs[])
{
	double t, r, th, phi;
	double dt, dr, dth, dphi;

	t = vars[0];
	r = vars[1];
	th = vars[2];
	phi = vars[3];

	dt = vars[4];
	dr = vars[5];
	dth = vars[6];
	dphi = vars[7];

	/* 1st order diff eqs */
	diffs[0] = dt;	 // dt
	diffs[1] = dr;	 // dr
	diffs[2] = dth;	 // dth
	diffs[3] = dphi; // dphi

	/* 2nd order diff eqs for r and theta; c.f. Eq (28) and Eq (29) */
	diffs[4] = -christ[0][0][0] * dt * dt - 2.0 * christ[0][0][1] * dt * dr - 2.0 * christ[0][0][2] * dt * dth - 2.0 * christ[0][0][3] * dt * dphi - christ[0][1][1] * dr * dr - 2.0 * christ[0][1][2] * dr * dth - 2.0 * christ[0][1][3] * dr * dphi - christ[0][2][2] * dth * dth - 2.0 * christ[0][2][3] * dth * dphi - christ[0][3][3] * dphi * dphi; // d2t
	diffs[5] = -christ[1][0][0] * dt * dt - 2.0 * christ[1][0][1] * dt * dr - 2.0 * christ[1][0][2] * dt * dth - 2.0 * christ[1][0][3] * dt * dphi - christ[1][1][1] * dr * dr - 2.0 * christ[1][1][2] * dr * dth - 2.0 * christ[1][1][3] * dr * dphi - christ[1][2][2] * dth * dth - 2.0 * christ[1][2][3] * dth * dphi - christ[1][3][3] * dphi * dphi; // d2r
	diffs[6] = -christ[2][0][0] * dt * dt - 2.0 * christ[2][0][1] * dt * dr - 2.0 * christ[2][0][2] * dt * dth - 2.0 * christ[2][0][3] * dt * dphi - christ[2][1][1] * dr * dr - 2.0 * christ[2][1][2] * dr * dth - 2.0 * christ[2][1][3] * dr * dphi - christ[2][2][2] * dth * dth - 2.0 * christ[2][2][3] * dth * dphi - christ[2][3][3] * dphi * dphi; // d2th
	diffs[7] = -christ[3][0][0] * dt * dt - 2.0 * christ[3][0][1] * dt * dr - 2.0 * christ[3][0][2] * dt * dth - 2.0 * christ[3][0][3] * dt * dphi - christ[3][1][1] * dr * dr - 2.0 * christ[3][1][2] * dr * dth - 2.0 * christ[3][1][3] * dr * dphi - christ[3][2][2] * dth * dth - 2.0 * christ[3][2][3] * dth * dphi - christ[3][3][3] * dphi * dphi; // d2phi
}

void diffeqs(double vars[], double diffs[])
{
	double t, r, th, phi;
	double dt, dr, dth, dphi;
	double christ[4][4][4];

	t = vars[0];
	r = vars[1];
	th = vars[2];
	phi = vars[3];

	dt = vars[4];
	dr = vars[5];
	dth = vars[6];
	dphi = vars[7];

	christoffel(spin, defpar, r, th, christ);

	/* 1st order diff eqs */
	diffs[0] = dt;	 // dt
	diffs[1] = dr;	 // dr
	diffs[2] = dth;	 // dth
	diffs[3] = dphi; // dphi

	/* 2nd order diff eqs for r and theta; c.f. Eq (28) and Eq (29) */
	diffs[4] = -christ[0][0][0] * dt * dt - 2.0 * christ[0][0][1] * dt * dr - 2.0 * christ[0][0][2] * dt * dth - 2.0 * christ[0][0][3] * dt * dphi - christ[0][1][1] * dr * dr - 2.0 * christ[0][1][2] * dr * dth - 2.0 * christ[0][1][3] * dr * dphi - christ[0][2][2] * dth * dth - 2.0 * christ[0][2][3] * dth * dphi - christ[0][3][3] * dphi * dphi; // d2t
	diffs[5] = -christ[1][0][0] * dt * dt - 2.0 * christ[1][0][1] * dt * dr - 2.0 * christ[1][0][2] * dt * dth - 2.0 * christ[1][0][3] * dt * dphi - christ[1][1][1] * dr * dr - 2.0 * christ[1][1][2] * dr * dth - 2.0 * christ[1][1][3] * dr * dphi - christ[1][2][2] * dth * dth - 2.0 * christ[1][2][3] * dth * dphi - christ[1][3][3] * dphi * dphi; // d2r
	diffs[6] = -christ[2][0][0] * dt * dt - 2.0 * christ[2][0][1] * dt * dr - 2.0 * christ[2][0][2] * dt * dth - 2.0 * christ[2][0][3] * dt * dphi - christ[2][1][1] * dr * dr - 2.0 * christ[2][1][2] * dr * dth - 2.0 * christ[2][1][3] * dr * dphi - christ[2][2][2] * dth * dth - 2.0 * christ[2][2][3] * dth * dphi - christ[2][3][3] * dphi * dphi; // d2th
	diffs[7] = -christ[3][0][0] * dt * dt - 2.0 * christ[3][0][1] * dt * dr - 2.0 * christ[3][0][2] * dt * dth - 2.0 * christ[3][0][3] * dt * dphi - christ[3][1][1] * dr * dr - 2.0 * christ[3][1][2] * dr * dth - 2.0 * christ[3][1][3] * dr * dphi - christ[3][2][2] * dth * dth - 2.0 * christ[3][2][3] * dth * dphi - christ[3][3][3] * dphi * dphi; // d2phi
}


int raytrace(double inc, double spin, double defpar, double dscr, double xscr, double yscr, double rin, double rout, double traced[])
{
    double xscr2, yscr2;
    double t, r, th, phi, tau, rau, thau, phiau;
    double t0, r0, th0, phi0;
    double kt, kr, kth, kphi;
    double kt0, kr0, kth0, kphi0;
    double s0, r02, s02;
    double fact1, fact2, fact3, B, C, omega;
    double h;

	double dd, ss, ssss, horizon;

    double height;

    double cosem;
    double b;

    double g[4][4]; // metric tensor
    double diffs[8], vars[8], vars_temp[8], vars_4th[8], vars_5th[8], k1[8], k2[8], k3[8], k4[8], k5[8], k6[8];
	double xem[4];
    double rgap, rmid, thmid;
    double gfactor;
    double err, errmin, errmax;
	double err_offset;
    double cross_tol;

    int errcheck, crosscheck = 0, acccheck = 0, blockcheck = 0;
    int i;
	int stop_integration = 0;

    /* ----- Set computational parameters ----- */
	err_offset = 1.0e3;
    errmin = 1.0e-9; // error bounds for RK45
    errmax = 1.0e-7;
    cross_tol = 1.0e-8; // sought accuracy at disk crossing

    // Set disk height
    // if (rdisk >= isco)
    //     height = 3.0 * Mdl * (1.0 - sqrt(isco / rdisk)) / eta;
    // else
    //     height = 0.0;
	height = 0.0;

    h = -1.0; // initial step size, negative sign for backwards integration

	double spin2 = spin*spin;

    /* ----- Compute photon initial conditions ----- */

    // Variables used for convenience
    xscr2 = xscr * xscr;
    yscr2 = yscr * yscr;
    fact1 = yscr * sin(inc) + dscr * cos(inc); // c.f. numerator of Eq (31)
    fact2 = dscr * sin(inc) - yscr * cos(inc); // c.f. denominator of Eq(32)
    r02 = xscr2 + yscr2 + dscr * dscr;         // c.f. Eq (30)

    // Initial t, r, theta, and phi; c.f. Eqs (30)-(32)
    t0 = 0.0;
    r0 = sqrt(r02); 
    th0 = acos(fact1 / r0);
    phi0 = atan2(xscr, fact2); // atan2 used to get principal value

    // Initial r, theta, and phi momentum; c.f. Eqs (33)-(35)
    kr0 = dscr / r0;
    kth0 = -(cos(inc) - dscr * fact1 / r02) / sqrt(r02 - fact1 * fact1);
    kphi0 = -xscr * sin(inc) / (xscr2 + fact2 * fact2);

    metric(r0, th0, g); // compute initial metric tensor
    fact3 = sqrt(g[0][3] * g[0][3] * kphi0 * kphi0 - g[0][0] * (g[1][1] * kr0 * kr0 + g[2][2] * kth0 * kth0 + g[3][3] * kphi0 * kphi0));
    B = 2.0 * g[0][1] * kr0 + 2 * g[0][3] * kphi0;                                                         // B
    C = g[1][1] * kr0 * kr0 + 2 * g[1][3] * kr0 * kphi0 + g[2][2] * kth0 * kth0 + g[3][3] * kphi0 * kphi0; // C

    // Initial t momentum
    kt0 = (B - sqrt(B * B - 4.0 * C * g[0][0])) / (2.0 * g[0][0]);

    // Impact parameter, b=L/E
    b = -(g[3][3] * kphi0 + g[0][3] * kt0) / (g[0][0] * kt0 + g[0][3] * kphi0);

    // Scale by Energy
    kr0 /= fact3;
    kth0 /= fact3;

    /* ----- RK45 ----- */

    /* set initial values */
    t = t0;
    r = r0;
    th = th0;
    phi = phi0;

    kt = kt0;
    kr = kr0;
    kth = kth0;
    kphi = kphi0;

    s0 = sin(th0);
    s02 = s0 * s0;

    omega = r02 * s02 * kphi0 / kt0; /* disk angular velocity */

    do
    {

        vars[0] = t;
        vars[1] = r;
        vars[2] = th;
        vars[3] = phi;
        vars[4] = kt;
        vars[5] = kr;
        vars[6] = kth;
        vars[7] = kphi;

        do
        {

            errcheck = 0;

            /* ----- compute RK1 ----- */

            diffeqs(vars, diffs);
            for (i = 0; i <= 7; i++)
            {
                k1[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a1 * k1[i];
            }

            /* ----- compute RK2 ----- */

            diffeqs(vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k2[i] = h * diffs[i];
                vars_temp[i] = vars[i] + b1_rk * k1[i] + b2_rk * k2[i];
            }

            /* ----- compute RK3 ----- */

            diffeqs(vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k3[i] = h * diffs[i];
                vars_temp[i] = vars[i] + c1_rk * k1[i] + c2_rk * k2[i] + c3 * k3[i];
            }

            /* ----- compute RK4 ----- */

            diffeqs(vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k4[i] = h * diffs[i];
                vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i] + d4 * k4[i];
            }

            /* ----- compute RK5 ----- */

            diffeqs(vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k5[i] = h * diffs[i];
                vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i];
            }

            /* ----- compute RK6 ----- */

            diffeqs(vars_temp, diffs);
            for (i = 0; i <= 7; i++)
                k6[i] = h * diffs[i];

            /* ----- local error ----- */

            for (i = 0; i <= 7; i++)
            {
                vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i] + f4 * k4[i] + f5 * k5[i];
                vars_5th[i] = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i] + g4 * k4[i] + g5 * k5[i] + g6 * k6[i];

                if (i == 0 || i == 4)
                {
                      ;
                } // pass
                else
                {
                err = fabs((vars_4th[i] - vars_5th[i]) / max(vars_4th[i], vars[i]));

                if (err > errmax) /* accuracy not achieved and photon hasn't crossed disk */
                    errcheck = 1;
                else if (err < errmin && errcheck != 1) /* accuracy better than wanted, but photon hasn't crossed disk */
                    errcheck = -1;
            	}
			}

            if (errcheck == 1) /* accuracy not achieved, lower step size */
                h /= 2.0;
			else if (errcheck == -1) /* accuracy better than wanted, increase step size */
				h *= 2.0;

        } while (errcheck == 1);

        /* ----- cross disk/horizon check ----- */

        tau = t;
        rau = r;
        thau = th;
        phiau = phi;
        t = vars_4th[0];
        r = vars_4th[1];
        th = vars_4th[2];
        phi = vars_4th[3];
		kt = vars_4th[4];
		kr = vars_4th[5];
		kth = vars_4th[6];
		kphi = vars_4th[7];

        if (cos(th) < 0.0) // check if photon has crossed disk
        {
			intersection(rau, thau, phiau, r, th, phi, xem);
			
            /* check if accuracy achieved; if so, move on to setting final values */
            if (sqrt(r * r + rau * rau - 2.0 * r * rau * (cos(th) * cos(thau) + sin(th) * sin(thau) * cos(phi - phiau))) <= cross_tol)
            {
                acccheck = 1;
                // printf("Crossed disk near desired radius and within error tolerance. Exit.\n");
            }
            /* otherwise, photon has crossed disk, but has not achieved accuracy; continue to zoom in on disk */
            else
            {
                crosscheck = 1;
                // printf("Crossed disk near desired radius but not within error tolerance. Continue.\n");
            }
            
            if (acccheck == 1) /* set final values */
            {
                /* calculate average/midpoint values */
                rmid = 0.5 * (r + rau);
                thmid = 0.5 * (th + thau);

                printf("%Le %Le\n", rmid, thmid);

                if (rmid * sin(thmid) >= rin - 0.001 && rmid * sin(thmid) < 1.05 * dscr)
                {
                    stop_integration = 1; /* the photon hits the disk */
                    break;
                }
                else
                {
                    stop_integration = 2; /* the photon misses the disk or other error */
                    break;
                }
            }
            else if (crosscheck == 1) /* did not achieve accuracy; go back a step, and decrease step size */
            {
                r = rau;
                th = thau;
                phi = phiau;
                h /= 2.0;
            }
        }
        else if (r <= 1. + sqrt(1. - spin2) + 0.001) // check if photon has crossed horizon
        {
            stop_integration = 2;
            // printf("Photon crossed horizon\n");
            break;
        }
        else if (h > -1.0e-20)
        {
            stop_integration = 3; /* photon is stuck */
            // printf("Photon is stuck\n");
            break;
        }
        else /* not done, take a step */
        {
			;
        }
    } while (stop_integration == 0);

    /* ----- Calculate redshift, cosem, and return values ----- */

    if (stop_integration == 1) /* photon hit disk, no issues */
    {
		intersection(rau, thau, phiau, r, th, phi, xem);
        redshift(spin, defpar, xem[1], th, omega, gfactor);
    }
    else /* photon crossed horizon, missed disk, or other issue */
    {
        rmid = 0.0;
        gfactor = 0.0;
    }

    traced[0] = rmid;
    traced[1] = gfactor;

	return stop_integration;
}
