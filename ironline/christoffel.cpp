void christoffel(double r, double th, double christ[4][4][4])
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