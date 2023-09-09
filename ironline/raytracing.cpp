int raytrace(double errmin, double errmax, double xscr, double yscr, double traced[4])
{
    double dscr, xscr2, yscr2;
    double t, r, th, phi, tau, rau, thau, phiau;
    double t0, r0, th0, phi0;
    double kt, kr, kth, kphi, ktau, krau, kthau, kphiau;
    double kt0, kr0, kth0, kphi0;
    double s0, r02, s02;
    double fact1, fact2, fact3, const1;
    double h;
    double xem[4];

    double height;

    double g[4][4]; // metric tensor
    double diffs[8], vars[8], vars_temp[8], vars_4th[8], vars_5th[8], k1[8], k2[8], k3[8], k4[8], k5[8], k6[8];
    double rmid, thmid;
    double gfactor, limbdark;
    double err, errtol;
    double cross_tol;

    int errcheck;
    int stop_integration = 0;
    int i;

    /* ----- Set computational parameters ----- */
    dscr = 1.0e+6;   // distance of the observer, effectively at infinity
    errtol = 1.0e-8; // error tolerance for adaptive step size
    cross_tol = 1.0e-8; // sought accuracy at disk crossing

    height = 0.0;

    h = -1.0; // initial step size, negative sign for backwards integration

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
    kt0 = -sqrt(kr0 * kr0 + r02 * kth0 * kth0 + r02 * sin(th0) * sin(th0) * kphi0 * kphi0);

    // // Impact parameter, b=L/E
    // b = -(g[3][3] * kphi0 + g[0][3] * kt0) / (g[0][0] * kt0 + g[0][3] * kphi0);

    // // Scale by Energy
    // kr0 /= fact3;
    // kth0 /= fact3;

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

    const1 = r02 * s02 * kphi0 / kt0; /* disk angular velocity */

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
                }

                // err = fabs((vars_4th[i] - vars_5th[i]) / max(vars_4th[i], vars[i]));

                if (err > errmax) /* accuracy not achieved and photon hasn't crossed disk */
                    errcheck = 1;
                else if (err < errmin) /* accuracy better than wanted, but photon hasn't crossed disk */
                    errcheck = -1;
            
            }

            if (errcheck == 1) /* accuracy not achieved, lower step size */
                h /= 2.0;
            else if (errcheck == -1) /* accuracy better than wanted, but photon hasn't crossed disk, increase step size */
                h *= 2.0;

        } while (errcheck == 1);

        /* ----- cross disk/horizon check ----- */

        tau = t;
        rau = r;
        thau = th;
        phiau = phi;
        ktau = kt;
        krau = kr;
        kthau = kth;
        kphiau = kphi;

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
            rmid = (rau + r) / 2.0;
            thmid = (thau + th) / 2.0;

            intersection(rau, thau, phiau, r, th, phi, xem);

            if (xem[1] > isco && xem[1] < isco+250.0) {
                stop_integration = 1;   /* the photon hits the disk */
            } else {
                stop_integration = 2;   /* the photon misses the disk */
            }

            // printf("%Le %Le\n", rmid, thmid);

            // if (rmid * sin(thmid) >= isco - 0.001 && rmid * sin(thmid) < 1.05 * dscr)
            // {
            //     stop_integration = 1; /* the photon hits the disk */
            //     intersection(rau, thau, phiau, r, th, phi, xem);
            // }
            // else
            // {
            //     stop_integration = 2; /* the photon misses the disk or other error */
            // }
        }
        else if (r <= 1. + sqrt(1. - spin2) + 0.001) // check if photon has crossed horizon
        {
            stop_integration = 2;
            // printf("Photon crossed horizon\n");
        }
        else if (h > -1.0e-20)
        {
            stop_integration = 3; /* photon is stuck */
            // printf("Photon is stuck\n");
        }
        else if (r >= 1.05*dobs)
        {
            stop_integration = 4; /* photon is too far away */
            // printf("Photon is too far away\n");
        }
        else /* not done, take a step */
        {
            ;
        }
    } while (stop_integration == 0);

    /* ----- Calculate redshift and return values ----- */

    if (stop_integration == 1) /* photon hit disk, no issues */
    {
        redshift(rmid, Pi/2, const1, gfactor);
    }
    else /* photon crossed horizon, missed disk, or other issue */
    {
        rmid = 0.0;
        gfactor = 0.0;
    }

    traced[0] = rmid;
    traced[1] = gfactor;
    traced[2] = xem[1];

    return stop_integration;
}


int raytrace_RKN(double errtol, double xobs, double yobs, double traced[4])
{
	double robs, pobs;
	double robs_i, robs_f, rstep, rstep2, pstep;
	double xin, xout;
	double t0, r0, th0, phi0;
	double kt0, kr0, kth0, kphi0;
	double r02, s0, s02;
	double fact1, fact2, fact3, fact4;
	double t, r, th, phi;
	double kt, kr, kth, kphi;
	double tau, rau, thau, phiau;
    double ktau, krau, kthau, kphiau;
	double xmap, ymap;
	double const0, const1;
	double h, hnext;
	double v1, v2;
	double check;
	double xdist;
	double dd, ss, ssss, horizon;
	double deltax, deltay, dxdy;
	double pp, qq;
	double fr;
	double Upsilon;

	double isco;
	double Gamma[4][4][4];
	double v[4], p[4];
	double u[4];
	double RK1[4], RK2[4], RK3[4], RK4[4];
	double verr[4], vtol[4], err[4];
	double xmid, xem[4];
	double gfactor;
	double limbdark;

    int i, n1, n2, n3;
    int stop_integration;
    
    r02 = xobs * xobs + yobs * yobs + dobs * dobs;
    r0 = sqrt(r02);

    fact1 = dobs * sin(inc) - yobs * cos(inc);
    fact2 = fact1 * fact1;
    fact3 = xobs * xobs + fact2;
    fact4 = sqrt(fact3);

    t0 = 0.0;
    r0 = r0;
    th0 = acos((yobs * sin(inc) + dobs * cos(inc)) / r0);
    phi0 = atan(xobs / fact1);

    s0 = sin(th0);
    s02 = s0 * s0;

    kr0 = -dobs / r0;
    kth0 = (cos(inc) - dobs * (yobs * sin(inc) + dobs * cos(inc)) / r02) / fact4;
    kphi0 = xobs * sin(inc) / fact3;
    kt0 = sqrt(kr0 * kr0 + r02 * kth0 * kth0 + r02 * s02 * kphi0 * kphi0);

    /* ----- solve geodesic equations
        fourth-order runge-kutta-nystrom method
        see E. Lund et al., 2009 JINST 4 P04001 ----- */

    t = t0;
    r = r0;
    th = th0;
    phi = phi0;

    kt = kt0;
    kr = kr0;
    kth = kth0;
    kphi = kphi0;

    const0 = kt0;
    const1 = r02 * s02 * kphi0 / kt0;

    stop_integration = 0;

    hnext = 100.0;

    do
    {

        v[0] = t;
        v[1] = r;
        v[2] = th;
        v[3] = phi;

        p[0] = kt;
        p[1] = kr;
        p[2] = kth;
        p[3] = kphi;

        do
        {

            h = hnext;

            /* ----- compute RK1 ----- */

            v1 = r;
            v2 = th;

            // christoffel(v1, v2, Gamma);
            Christoffel_jiale(v1, v2, Gamma);

            for (i = 0; i <= 3; i++)
                u[i] = p[i];

            for (n1 = 0; n1 <= 3; n1++)
            {
                for (n2 = 0; n2 <= 3; n2++)
                {
                    for (n3 = 0; n3 <= 3; n3++)
                    {
                        if (n2 == 0 && n3 == 0)
                        {
                            RK1[n1] = -Gamma[n1][0][0] * u[0] * u[0];
                        }
                        else
                        {
                            RK1[n1] -= Gamma[n1][n2][n3] * u[n2] * u[n3];
                        }
                    }
                }
            }

            /* ----- compute RK2 ----- */

            v1 = r + h * kr / 2 + h * h * RK1[1] / 16;
            v2 = th + h * kth / 2 + h * h * RK1[2] / 16;
            // v1 = r + h * kr / 2 + h * h * RK1[1] / 8;
            // v2 = th + h * kth / 2 + h * h * RK1[2] / 8;

            // christoffel(v1, v2, Gamma);
            Christoffel_jiale(v1, v2, Gamma);

            for (i = 0; i <= 3; i++)
                u[i] = p[i] + h * RK1[i] / 4;
                // u[i] = p[i] + h * RK1[i] / 2;

            for (n1 = 0; n1 <= 3; n1++)
            {
                for (n2 = 0; n2 <= 3; n2++)
                {
                    for (n3 = 0; n3 <= 3; n3++)
                    {
                        if (n2 == 0 && n3 == 0)
                        {
                            RK2[n1] = -Gamma[n1][0][0] * u[0] * u[0];
                        }
                        else
                        {
                            RK2[n1] -= Gamma[n1][n2][n3] * u[n2] * u[n3];
                        }
                    }
                }
            }

            /* ----- compute RK3 ----- */

            for (i = 0; i <= 3; i++)
                u[i] = p[i] + h * RK2[i] / 4;
                // u[i] = p[i] + h * RK1[i] / 2;

            for (n1 = 0; n1 <= 3; n1++)
            {
                for (n2 = 0; n2 <= 3; n2++)
                {
                    for (n3 = 0; n3 <= 3; n3++)
                    {
                        if (n2 == 0 && n3 == 0)
                        {
                            RK3[n1] = -Gamma[n1][0][0] * u[0] * u[0];
                        }
                        else
                        {
                            RK3[n1] -= Gamma[n1][n2][n3] * u[n2] * u[n3];
                        }
                    }
                }
            }

            /* ----- compute RK4 ----- */

            v1 = r + h * kr + h * h * RK3[1] / 4;
            v2 = th + h * kth + h * h * RK3[2] / 4;
            // v1 = r + h * kr + h * h * RK3[1] / 2;
            // v2 = th + h * kth + h * h * RK3[2] / 2;

            // christoffel(v1, v2, Gamma);
            Christoffel_jiale(v1, v2, Gamma);

            for (i = 0; i <= 3; i++)
                u[i] = p[i] + h * RK3[i] / 2;
                // u[i] = p[i] + h * RK3[i];

            for (n1 = 0; n1 <= 3; n1++)
            {
                for (n2 = 0; n2 <= 3; n2++)
                {
                    for (n3 = 0; n3 <= 3; n3++)
                    {
                        if (n2 == 0 && n3 == 0)
                        {
                            RK4[n1] = -Gamma[n1][0][0] * u[0] * u[0];
                        }
                        else
                        {
                            RK4[n1] -= Gamma[n1][n2][n3] * u[n2] * u[n3];
                        }
                    }
                }
            }

            /* ----- local error ----- */

            for (i = 0; i <= 3; i++)
            {
                verr[i] = 0.5 * h * h * (RK1[i] - RK2[i] - RK3[i] + RK4[i]);
                verr[i] *= verr[i];
                vtol[i] = errtol + fabs(v[i]) * errtol;
                vtol[i] *= vtol[i];
                err[i] = 0.25 * verr[i] / vtol[i];
            }

            check = sqrt(err[0] + err[1] + err[2] + err[3]);

            /* ----- next step ----- */

            hnext = h * pow(1 / check, 0.25);

            /* ----- limitation criterion ----- */

            if (hnext < h / 4)
                hnext = h / 4;
            if (hnext > 4 * h)
                hnext = 4 * h;

            if (hnext == h)
                hnext = 0.9 * h;

        } while (check > 1);

        /* ----- solutions to the fourth-order RKN method ----- */

        tau = t;
        rau = r;
        thau = th;
        phiau = phi;

        ktau = kt;
        krau = kr;
        kthau = kth;
        kphiau = kphi;

        t += h * kt + (RK1[0] + RK2[0] + RK3[0]) * h * h / 12;
        r += h * kr + (RK1[1] + RK2[1] + RK3[1]) * h * h / 12;
        th += h * kth + (RK1[2] + RK2[2] + RK3[2]) * h * h / 12;
        phi += h * kphi + (RK1[3] + RK2[3] + RK3[3]) * h * h / 12;

        kt += (RK1[0] + 2 * RK2[0] + 2 * RK3[0] + RK4[0]) * h / 12;
        kr += (RK1[1] + 2 * RK2[1] + 2 * RK3[1] + RK4[1]) * h / 12;
        kth += (RK1[2] + 2 * RK2[2] + 2 * RK3[2] + RK4[2]) * h / 12;
        kphi += (RK1[3] + 2 * RK2[3] + 2 * RK3[3] + RK4[3]) * h / 12;

        // t += h * kt + (RK1[0] + RK2[0] + RK3[0]) * h * h / 6;
        // r += h * kr + (RK1[1] + RK2[1] + RK3[1]) * h * h / 6;
        // th += h * kth + (RK1[2] + RK2[2] + RK3[2]) * h * h / 6;
        // phi += h * kphi + (RK1[3] + RK2[3] + RK3[3]) * h * h / 6;

        // kt += (RK1[0] + 2 * RK2[0] + 2 * RK3[0] + RK4[0]) * h / 6;
        // kr += (RK1[1] + 2 * RK2[1] + 2 * RK3[1] + RK4[1]) * h / 6;
        // kth += (RK1[2] + 2 * RK2[2] + 2 * RK3[2] + RK4[2]) * h / 6;
        // kphi += (RK1[3] + 2 * RK2[3] + 2 * RK3[3] + RK4[3]) * h / 6;

        if (cos(th) < 0.0)
        {
            intersection(rau, thau, phiau, r, th, phi, xem);

            if (xem[1] > isco - 0.001 && xem[1] < isco + 250)
            {
                stop_integration = 1; /* the photon hits the disk */
            }
            else
            {
                stop_integration = 2; /* the photon misses the disk */
            }
        }

        // if (r <= 1. + sqrt(1. - spin2) + 0.001)
        if (1 + sqrt(1 - spin2) < 0.001)
            stop_integration = 4; /* the photon crosses the horizon */

        if (r < 1)
            stop_integration = 5; /* the photon hits the singularity */

        if (r > 1.05 * dobs)
            stop_integration = 8; /* the photon escapes to infinity */

    } while (stop_integration == 0);

    /* ----- Calculate redshift and return values ----- */
    if (stop_integration == 1) /* photon hit disk, no issues */
    {
        redshift(xem[1], Pi / 2, const1, gfactor);
    }
    else /* photon crossed horizon, missed disk, or other issue */
    {
        xem[1] = 0.0;
        gfactor = 0.0;
        xmid = 0.0;
    }

    traced[0] = xmid;
    traced[1] = gfactor;
    traced[2] = xem[1];

    return stop_integration;
}