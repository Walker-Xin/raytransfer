#include "def.h"

int main(int argc, char *argv[])
{
	double cosinc;
	double rdisk_i, rdisk_f, rdisk[imax + 2];
	double traced[4];
	double gerrtol, rerrtol, pdiff;

	int mu_len, rdisk_len;

	double time_taken, iteration_time, expected_time;
	int progress_check, skip, skip_end;

	char filename_o[128];
	FILE *foutput;

	/* Set default computational values */
	spin = 0.9;
	Mdl = 0;
	defpar = 0;		  // default = 0
	gerrtol = 1.0e-6; // default = 1.0e-6
	rerrtol = 1.0e-7; // default = 1.0e-7
	pdiff = 1.0e-4;	  // default = 1.0e-4
	progress_check = 1;
	skip = 0;	   // default = 0
	skip_end = imax - 1; // default = imax - 1

	// Set computation parameters from user input if provided
	if (argc > 1)
	{
		spin = atof(argv[1]); // spin parameter
		inc = atof(argv[2]);  // inclination angle in degrees
		defpar = atof(argv[3]);
		gerrtol = atof(argv[4]);		// error tolerance for RK45
		rerrtol = atof(argv[5]);		// error tolerance for redshift factor
		progress_check = atoi(argv[6]); // progress check interval, also used to determine whether to print progress
		printf("Using user input parameters. spin=%.5e, inc=%.5e, deformation=%.5e, gerrtol=%.5e, rerrtol=%.5e\n", double(spin), double(inc), double(defpar), double(gerrtol), double(rerrtol));
	}
	else
	{
		printf("Using user input parameters. spin=%.5e, inc=%.5e, deformation=%.5e, gerrtol=%.5e, rerrtol=%.5e\n", double(spin), double(inc), double(defpar), double(gerrtol), double(rerrtol));
	}

	spin2 = spin * spin;

	/* Preset g_star values */
	double g_star[40] = {0.002, 0.02753846, 0.05307692, 0.07861538, 0.10415385, 0.12969231, 0.15523077, 0.18076923, 0.20630769, 0.23184615, 0.25738462, 0.28292308, 0.30846154, 0.334, 0.35953846, 0.38507692, 0.41061538, 0.43615385, 0.46169231, 0.48723077, 0.51276923, 0.53830769, 0.56384615, 0.58938462, 0.61492308, 0.64046154, 0.666, 0.69153846, 0.71707692, 0.74261538, 0.76815385, 0.79369231, 0.81923077, 0.84476923, 0.87030769, 0.89584615, 0.92138462, 0.94692308, 0.97246154, 0.998};

	/* Loop over inclination angles when running on cluster */
	// double mu0[] = {0.0349447653, 0.09718278, 0.15948, 0.2165542, 0.270481, 0.3221819, 0.3721757, 0.420793, 0.4682622, 0.5147499, 0.5603828, 0.6052601, 0.6494616, 0.6930526, 0.7360878, 0.7786132, 0.8206683, 0.8622873, 0.9035001, 0.9443328, 0.9848086238, 0.9986296296};
	double mu0[] = {cos(inc)};
	mu_len = sizeof(mu0) / sizeof(mu0[0]);

	/* Set inner radius of the disk */
	isco = find_isco();
	printf("isco = %.15Le\n", double(isco));

	/* Calculate radiative efficiency */
	eta = 1.0 - specific_energy(isco);
	printf("eta = %.15Le\n", double(eta));

	// Start timer
	clock_t start, end, mid;
	start = clock();
	mid = clock();

	printf("\nSIMULATION START\n");
	printf("----------------\n\n");

	for (int jj = 0; jj < mu_len; jj++)
	{
		cosinc = mu0[jj];
		inc = acos(cosinc); // inclination angle of the observer in rad
		printf("spin = %.5e, inc = %.5e (rad)/%.5e (deg) defpar = %.5e\n", double(spin), double(inc), double(inc * 180.0 / Pi), double(defpar));

		// /* Set inner radius of the disk */
		// isco = find_isco();
		// printf("isco = %.15Le\n", isco);

		// /* Calculate radiative efficiency */
		// eta = 1.0 - specific_energy(isco);
		// printf("eta = %.15Le\n", eta);

		/* Set inner/outer disk radii */
		rdisk_i = isco;
		rdisk_f = 1000.;

		/* Compute emission radii with gauleg */
		gauleg(rdisk_i, rdisk_f, rdisk); // a grid of 100 emission radii stored in rdisk; c.f. Section 3.4 in Public Release
		rdisk_len = sizeof(rdisk) / sizeof(rdisk[0]);

		/* Open output file */
		sprintf(filename_o, "photons/photons4trf_a%.05Le.i%.02Le.Mdl_%.02Le.dp_%.02Le.dat", double(spin), double(cosinc), double(Mdl), double(defpar));

		foutput = fopen(filename_o, "w");

		/* Assign photon position in the grid */
		for (int ii = skip; ii <= skip_end; ii++)
		{
			// !!NB!!
			// Redefining the variables at each loop is necessary for resetting the values
			double pscr, pstep, pscrcur, pscrhigh, pscrlow, rscr, rdiskcur, cosem, gcur, gstar, pdifft, gerrttol;
			double xscrcur, yscrcur, xscrplus, xscrminus, yscrplus, yscrminus, gplus, gminus;
			double gmax = 0.0, gmin = 10.0, pscrmax, pscrmin, rscrmax, rscrmin, cosemmax, cosemmin, rdiskmax, rdiskmin;

			if (progress_check != 0)
			{
				// Update progress for every progress_check robs
				if (ii % progress_check == 0)
				{
					// Calculate expected time usage
					iteration_time = double(clock() - mid) / double(CLOCKS_PER_SEC);
					mid = clock();
					expected_time = iteration_time * (rdisk_len - ii) / double(progress_check);
					printf("rdisk = %Le; expected time left: %Le minutes\n", double(rdisk[ii]), double(expected_time / 60.0));
				}
			}

			/* ------- Search over pscr to get quick estimate of gmin and gmax ------- */
			pstep = Pi / 5.0;
			pscr = 0.0;

			while (pscr < 2.0 * Pi)
			{
				rayprecise(rdisk[ii], rerrtol, pscr, traced);

				// printf("%Le %Le %Le\n", traced[2], gmin, gmax);

				if (traced[2] > gmax && traced[0] != 0.0) // set gmax if g is larger than current gmax
				{
					gmax = traced[2];
					pscrmax = pscr;
				}
				if (traced[2] < gmin && traced[0] != 0.0) // set gmin if g is smaller than current gmin
				{
					gmin = traced[2];
					pscrmin = pscr;
				}

				pscr += pstep;
			}

			/* ------- Search for gmax -------- */

			while (pstep > gerrtol)
			{
				rayprecise(rdisk[ii], rerrtol, pscrmax - pstep / 2.0, traced); // search at values of phi_screen lower than estimated

				// printf("%Le %Le %Le\n", traced[2], gmin, gmax);

				if (traced[2] > gmax)
				{
					gmax = traced[2];
					rdiskmax = traced[0];
					cosemmax = traced[1];
					rscrmax = traced[3];
					pscrmax -= pstep / 2.0;
				}

				// Can skip next rayprecise if gmax is updated?

				rayprecise(rdisk[ii], rerrtol, pscrmax + pstep / 2.0, traced); // search at values of phi_screen higher than estimated

				if (traced[2] > gmax)
				{
					gmax = traced[2];
					rdiskmax = traced[0];
					cosemmax = traced[1];
					rscrmax = traced[3];
					pscrmax += pstep / 2.0;
				}

				pstep /= 2.0;
			}

			/* ------- Search for gmin -------- */

			pstep = Pi / 5.0;

			while (pstep > gerrtol)
			{
				rayprecise(rdisk[ii], rerrtol, pscrmin - pstep / 2.0, traced); // search at values of phi_screen lower than estimated

				if (traced[2] < gmin)
				{
					gmin = traced[2];
					rdiskmin = traced[0];
					cosemmin = traced[1];
					rscrmin = traced[3];
					pscrmin -= pstep / 2.0;
				}

				// Can skip next rayprecise if gmin is updated?

				rayprecise(rdisk[ii], rerrtol, pscrmin + pstep / 2.0, traced); // search at values of phi_screen higher than estimated

				if (traced[2] < gmin)
				{
					gmin = traced[2];
					rdiskmin = traced[0];
					cosemmin = traced[1];
					rscrmin = traced[3];
					pscrmin += pstep / 2.0;
				}

				pstep /= 2.0;
			}

			if (progress_check != 0)
			{
				printf("%d gmin = %.15Le, gmax = %.15Le\n", ii, gmin, gmax);
			}

			//*_*_*_*_*_*_*_*		CALCULATING CONSTANTLY SPACED g* GRID		*_*_*_*_*_*_*//

			gerrttol = gerrtol;

			/*---------- Branch 1 ------------*/

			xyfromrphi(rscrmin, pscrmin, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", double(rdiskmin), double(gmin), double(xscr), double(yscr), double(cosemmin));

			if (progress_check != 0)
			{
				printf("%d B1:MIN %.6Le %.6Le %.6Le %.6Le %.6Le\n", ii, double(rdiskmin), double(gmin), double(xscr), double(yscr), double(cosemmin));
			}

			// set values for phi_screen variables
			pstep = fabs(pscrmax - pscrmin) / 39.0;
			pscrlow = pscrmin;
			pscrhigh = pscrmax;
			pscrcur = pscrmin + pstep; // searching from low to high

			for (int j = 0; j < 40; j++)
			{
				gcur = gmin + (gmax - gmin) * g_star[j];
				pstep = fabs(pscrmax - pscrmin) / 39.0;

				while (1) // loop until gcur is found
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur, traced);

					gstar = (traced[2] - gmin) / (gmax - gmin); // traced gstar value

					if (traced[2] == 0.0 || traced[0] == 0.0 || pscrcur >= pscrmax || pscrcur <= pscrmin) // Check that there's no error or going out of bounds of gmax and gmin
					{
						if (pscrcur >= pscrmax) // Lower phi_screen_current if past gmax value
						{
							pscrcur = pscrmax - pstep;
							pstep /= 2.0;
						}
						else if (pscrcur <= pscrmin) // Increase phi_screen_current if past gmin value
							pscrcur = pscrmin + pstep * 1.5;
						else // Step down for other errors
							pscrcur -= pstep * 0.5;
					}
					else if (fabs(gstar - g_star[j]) < gerrttol) // Found value, accuracy reached
					{
						// printf("Accuracy reached\n");
						break;
					}
					else if (traced[2] < gcur) // Value below gcur
						pscrlow = pscrcur;
					else if (traced[2] > gcur) // Value above gcur
						pscrhigh = pscrcur;

					// Set phi_screen_current to midpoint of high and low values
					// Possible truncation error here, especially for higher values of j
					// For example, low = 2.8624560870942282, high = 2.8624560870942286 gives cur = 2.8624560870942286 which is incorrect
					pscrcur = 0.5 * (pscrlow + pscrhigh);

					// Fix by breaking if cur = high or low after taking midpoint
					if (pscrcur == pscrhigh) // some truncation error has occured
					{
						printf("pscrcur = pscrhigh, j = %d\n", j);
						printf("pscrcur = %.16Le, pscrhigh = %.16Le\n, pscrlow = %.16Le", pscrcur, pscrhigh, pscrlow);
						printf("breaking\n");
						break;
					}
					else if (pscrcur == pscrlow)
					{
						printf("pscrcur = pscrlow, j = %d\n", j);
						printf("pscrcur = %.16Le, pscrhigh = %.16Le\n, pscrlow = %.16Le", pscrcur, pscrhigh, pscrlow);
						printf("breaking\n");
						break;
					}

					if (pscrhigh - pscrlow < 1.0e-10) // Raise error tolerance if can't find the correct value
						gerrttol *= 2.0;
				}

				// Set/Reset variables as necessary after finding gcur
				pscrlow = pscrcur;
				pscrhigh = pscrmax;

				gerrttol = gerrtol;

				rdiskcur = traced[0];
				gcur = traced[2];
				cosem = traced[1];
				rscr = traced[3];

				xyfromrphi(rscr, pscrcur, rdisk[ii]);
				xscrcur = xscr;
				yscrcur = yscr;

				pdifft = pdiff;

				// Find rays with phi_screen slightly lower and slightly higher than that for found gcur values
				do
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur - pdifft, traced);
					xyfromrphi(traced[3], pscrcur - pdifft, rdisk[ii]);
					xscrminus = xscr;
					yscrminus = yscr;
					gminus = traced[2];
					pdifft *= 2.0;
				} while (gminus == 0.0);

				pdifft = pdiff;

				do
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur + pdifft, traced);
					xyfromrphi(traced[3], pscrcur + pdifft, rdisk[ii]);
					xscrplus = xscr;
					yscrplus = yscr;
					gplus = traced[2];
					pdifft *= 2.0;
				} while (gplus == 0.0);

				fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le\n", double(rdiskcur), double(gcur), double(xscrcur), double(yscrcur), double(cosem), double(gminus), double(xscrminus), double(yscrminus), double(gplus), double(xscrplus), double(yscrplus));

				if (progress_check != 0)
				{
					printf("%d B1:%d %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le\r", ii, j + 1, double(rdiskcur), double(gcur), double(xscrcur), double(yscrcur), double(cosem), double(gminus), double(xscrminus), double(yscrminus), double(gplus), double(xscrplus), double(yscrplus));
				}
			}

			xyfromrphi(rscrmax, pscrmax, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", double(rdiskmax), double(gmax), double(xscr), double(yscr), double(cosemmax));

			if (progress_check != 0)
			{
				printf("%d B1:MAX %.6Le %.6Le %.6Le %.6Le %.6Le\n", ii, double(rdiskmax), double(gmax), double(xscr), double(yscr), double(cosemmax));
			}

			/*---------- Branch 2 ------------*/

			/* Note that for Branch 2 some things are done in reverse as the values are negative */

			xyfromrphi(rscrmin, pscrmin, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", double(rdiskmin), double(gmin), double(xscr), double(yscr), double(cosemmin));

			if (progress_check != 0)
			{
				printf("%d B2:MIN %.6Le %.6Le %.6Le %.6Le %.6Le\n", ii, double(rdiskmin), double(gmin), double(xscr), double(yscr), double(cosemmin));
			}

			// set values for phi_screen variables
			pstep = fabs(pscrmax - pscrmin) / 39.0;
			pscrlow = pscrmin;
			pscrhigh = pscrmax - 2.0 * Pi;
			pscrcur = pscrmin - pstep;

			for (int j = 0; j < 40; j++)
			{
				gcur = gmin + (gmax - gmin) * g_star[j];
				pstep = fabs(pscrmax - pscrmin) / 39.0;

				while (1) // loop until gcur is found
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur, traced);

					gstar = (traced[2] - gmin) / (gmax - gmin); // traced gstar value

					if (traced[2] == 0.0 || traced[0] == 0.0 || pscrcur <= pscrmax - 2.0 * Pi || pscrcur >= pscrmin) // Check that there's no error or going out of bounds of gmax and gmin
					{
						if (pscrcur <= pscrmax - 2.0 * Pi) // Raise phi_screen_current if past gmax value
						{
							pscrcur = pscrmax - 2.0 * Pi + pstep;
							pstep /= 2.0;
						}
						else if (pscrcur >= pscrmin) // Lower phi_screen_current if past gmin value
							pscrcur = pscrmin - pstep * 1.5;
						else // Step up for other errors
							pscrcur += pstep * 0.5;
					}
					else if (fabs(gstar - g_star[j]) < gerrttol) // Found value, accuracy reached
					{
						// printf("Accuracy reached\n");
						break;
					}
					else if (traced[2] < gcur) // Value below gcur
						pscrlow = pscrcur;
					else if (traced[2] > gcur) // Value above gcur
						pscrhigh = pscrcur;

					pscrcur = 0.5 * (pscrlow + pscrhigh); // Set phi_screen_current to midpoint of high and low values

					if (pscrcur == pscrhigh) // some truncation error has occured
					{
						printf("pscrcur = pscrhigh, j = %d\n", j);
						printf("pscrcur = %.16Le, pscrhigh = %.16Le\n, pscrlow = %.16Le", pscrcur, pscrhigh, pscrlow);
						printf("breaking\n");
						break;
					}
					else if (pscrcur == pscrlow)
					{
						printf("pscrcur = pscrlow, j = %d\n", j);
						printf("pscrcur = %.16Le, pscrhigh = %.16Le\n, pscrlow = %.16Le", pscrcur, pscrhigh, pscrlow);
						printf("breaking\n");
						break;
					}

					if (pscrlow - pscrhigh < 1.0e-10) // Raise error tolerance if can't find the correct value
						gerrttol *= 2.0;
				}

				// Set/Reset variables as necessary after finding gcur

				pscrlow = pscrcur;
				pscrhigh = pscrmax - 2.0 * Pi;

				gerrttol = gerrtol;

				rdiskcur = traced[0];
				gcur = traced[2];
				cosem = traced[1];
				rscr = traced[3];

				xyfromrphi(rscr, pscrcur, rdisk[ii]);
				xscrcur = xscr;
				yscrcur = yscr;

				pdifft = pdiff;

				// Find rays with phi_screen slightly lower and slightly higher than that for found gcur values

				do
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur - pdifft, traced);
					xyfromrphi(traced[3], pscrcur - pdifft, rdisk[ii]);
					xscrminus = xscr;
					yscrminus = yscr;
					gminus = traced[2];
					pdifft *= 2.0;
				} while (gminus == 0.0);

				pdifft = pdiff;

				do
				{
					rayprecise(rdisk[ii], rerrtol, pscrcur + pdifft, traced);
					xyfromrphi(traced[3], pscrcur + pdifft, rdisk[ii]);
					xscrplus = xscr;
					yscrplus = yscr;
					gplus = traced[2];
					pdifft *= 2.0;
				} while (gplus == 0.0);

				fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le\n", double(rdiskcur), double(gcur), double(xscrcur), double(yscrcur), double(cosem), double(gminus), double(xscrminus), double(yscrminus), double(gplus), double(xscrplus), double(yscrplus));

				if (progress_check != 0)
				{
					// Use /r to overwrite the line
					printf("%d B2:%d %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le\r", ii, j + 1, double(rdiskcur), double(gcur), double(xscrcur), double(yscrcur), double(cosem), double(gminus), double(xscrminus), double(yscrminus), double(gplus), double(xscrplus), double(yscrplus));
				}
			}

			xyfromrphi(rscrmax, pscrmax, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", double(rdiskmax), double(gmax), double(xscr), double(yscr), double(cosemmax));

			if (progress_check != 0)
			{
				printf("%d B2:MAX %.6Le %.6Le %.6Le %.6Le %.6Le\n", ii, double(rdiskmax), double(gmax), double(xscr), double(yscr), double(cosemmax));
			}
		}

		printf("----------------\n");
		printf("SIMULATION END\n");

		// Calculate total time taken
		end = clock();
		time_taken = double(end - start) / double(CLOCKS_PER_SEC);
		printf("Total time taken: %Le minutes\n", double(time_taken / 60.0));
		
		fclose(foutput);
	}

	return 0;
}
