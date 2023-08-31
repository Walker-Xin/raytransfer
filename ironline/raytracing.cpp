// TODO: change RK45 to RKN method

int raytrace(double xscr, double yscr, double traced[4])
{
    double dscr, xscr2, yscr2;
    double t, r, th, phi, tau, rau, thau, phiau;
    double t0, r0, th0, phi0;
    double kt, kr, kth, kphi, ktau, krau, kthau, kphiau;
    double kt0, kr0, kth0, kphi0;
    double s0, r02, s02;
    double fact1, fact2, fact3, B, C, omega;
    double h;
    double isco = 4.2330000000134405;
    double spin = 0.5;
    double spin2 = spin * spin;
    double xem[4];

    double height;

    double cosem;
    double b;

    double g[4][4]; // metric tensor
    double diffs[8], vars[8], vars_temp[8], vars_4th[8], vars_5th[8], k1[8], k2[8], k3[8], k4[8], k5[8], k6[8];
    double rgap, rmid, thmid;
    double gfactor, limbdark;
    double err, errmin, errmax, errtol;
    double cross_tol;

    int errcheck, crosscheck = 0, acccheck = 0, blockcheck = 0;
    int stop_integration = 0;
    int i;

    /* ----- Set computational parameters ----- */
    dscr = 1.0e+6;   // distance of the observer, effectively at infinity
    errtol = 1.0e-8; // error tolerance for adaptive step size
    errmin = 1.0e-11; // error bounds for RK45
    errmax = 1.0e-9;
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

            diffeqs(spin, defpar, vars, diffs);
            for (i = 0; i <= 7; i++)
            {
                k1[i] = h * diffs[i];
                vars_temp[i] = vars[i] + a1 * k1[i];
            }

            /* ----- compute RK2 ----- */

            diffeqs(spin, defpar, vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k2[i] = h * diffs[i];
                vars_temp[i] = vars[i] + b1_rk * k1[i] + b2_rk * k2[i];
            }

            /* ----- compute RK3 ----- */

            diffeqs(spin, defpar, vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k3[i] = h * diffs[i];
                vars_temp[i] = vars[i] + c1_rk * k1[i] + c2_rk * k2[i] + c3 * k3[i];
            }

            /* ----- compute RK4 ----- */

            diffeqs(spin, defpar, vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k4[i] = h * diffs[i];
                vars_temp[i] = vars[i] + d1 * k1[i] + d2 * k2[i] + d3 * k3[i] + d4 * k4[i];
            }

            /* ----- compute RK5 ----- */

            diffeqs(spin, defpar, vars_temp, diffs);
            for (i = 0; i <= 7; i++)
            {
                k5[i] = h * diffs[i];
                vars_temp[i] = vars[i] + e1 * k1[i] + e2 * k2[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i];
            }

            /* ----- compute RK6 ----- */

            diffeqs(spin, defpar, vars_temp, diffs);
            for (i = 0; i <= 7; i++)
                k6[i] = h * diffs[i];

            /* ----- local error ----- */

            for (i = 0; i <= 7; i++)
            {
                vars_4th[i] = vars[i] + f1 * k1[i] + f2 * k2[i] + f3 * k3[i] + f4 * k4[i] + f5 * k5[i];
                vars_5th[i] = vars[i] + g1 * k1[i] + g2 * k2[i] + g3 * k3[i] + g4 * k4[i] + g5 * k5[i] + g6 * k6[i];

                // if (i == 0 || i == 4)
                // {
                //       ;
                // } // pass
                // else
                // {
                // err = fabs((vars_4th[i] - vars_5th[i]) / max(vars_4th[i], vars[i]));
                // }

                err = fabs((vars_4th[i] - vars_5th[i]) / max(vars_4th[i], vars[i]));

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
            rmid = r;
            thmid = th;

            // printf("%Le %Le\n", rmid, thmid);

            if (rmid * sin(thmid) >= isco - 0.001 && rmid * sin(thmid) < 1.05 * dscr)
            {
                stop_integration = 1; /* the photon hits the disk */
                intersection(rau, thau, phiau, r, th, phi, xem);
                
                // // check xem[1] is the same as rmid
                // if (fabs(xem[1] - rmid) < 10e-8)
                // {
                //     break;
                // }
                // else
                // {
                //     // raise error
                //     printf("Error at r = %Le, th = %Le, phi = %Le\n", r, th, phi);
                //     exit(1);
                // }
            }
            else
            {
                stop_integration = 2; /* the photon misses the disk or other error */
                break;
            }
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
        else /* not done, take a step */
        {
            ;
        }
    } while (stop_integration == 0);

    /* ----- Calculate redshift and return values ----- */

    if (stop_integration == 1) /* photon hit disk, no issues */
    {
        redshift(spin, spin2, defpar, defpar, rmid, kt0, omega, kphiau, gfactor, limbdark);
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
