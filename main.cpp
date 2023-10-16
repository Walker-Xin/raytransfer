#include "def.h"

int main(int argc, char *argv[])
{
	double cosinc;
	double rdisk_i, rdisk_f, rdisk[imax + 2];
	double pscr, pstep, pscrcur, pscrhigh, pscrlow, rscr, rdiskcur, cosem, gcur, gstar, pdifft, gerrttol;
	double xscrcur, yscrcur, xscrplus, xscrminus, yscrplus, yscrminus, gplus, gminus;
	double gmax = 0.0, gmin = 10.0, pscrmax, pscrmin, rscrmax, rscrmin, cosemmax, cosemmin, rdiskmax, rdiskmin;
	double traced[4];
	double gerrtol, rerrtol, pdiff;

	int mu_len;

	double time_taken, iteration_time, expected_time;
	int progress_check;

	char filename_o[128];
	FILE *foutput;

	/* Set default computational values */
	spin = 0.1;
	Mdl = 0.1;
	defpar = 0.005;
	gerrtol = 1.0e-2;
	rerrtol = 1.0e-2;
	pdiff = 1.0e-4;
	progress_check = 1;

	// Set computation parameters from user input if provided
	if (argc > 1)
	{
		spin = atof(argv[1]); // spin parameter
		Mdl = atof(argv[2]);  // accretion rate parameter - disk thickness
		defpar = atof(argv[3]);
		gerrtol = atof(argv[4]); // error tolerance for RK45
		rerrtol = atof(argv[5]);
		pdiff = atof(argv[6]);
		printf("Using user input parameters. spin=%f, accretion rate=%f, deformation=%f, error tolerance=%e\n, rerrtol=%e, pdiff=%e\n", spin, Mdl, defpar, gerrtol, rerrtol, pdiff);
	}
	else
	{
		printf("Using preset parameters. spin=%f, accretion rate=%f, deformation=%f, error tolerance=%e\n, rerrtol=%e, pdiff=%e\n", spin, Mdl, defpar, gerrtol, rerrtol, pdiff);
	}

	spin2 = spin * spin;

	/* Preset g_star values */
	double g_star[40] = {0.002, 0.02753846, 0.05307692, 0.07861538, 0.10415385, 0.12969231, 0.15523077, 0.18076923, 0.20630769, 0.23184615, 0.25738462, 0.28292308, 0.30846154, 0.334, 0.35953846, 0.38507692, 0.41061538, 0.43615385, 0.46169231, 0.48723077, 0.51276923, 0.53830769, 0.56384615, 0.58938462, 0.61492308, 0.64046154, 0.666, 0.69153846, 0.71707692, 0.74261538, 0.76815385, 0.79369231, 0.81923077, 0.84476923, 0.87030769, 0.89584615, 0.92138462, 0.94692308, 0.97246154, 0.998};

	/* Loop over inclination angles when running on cluster */
	double mu0[] = {0.0349447653};
	// double mu0[] = {0.0349447653, 0.09718278, 0.15948, 0.2165542, 0.270481, 0.3221819, 0.3721757, 0.420793, 0.4682622, 0.5147499, 0.5603828, 0.6052601, 0.6494616, 0.6930526, 0.7360878, 0.7786132, 0.8206683, 0.8622873, 0.9035001, 0.9443328, 0.9848086238, 0.9986296296};
	mu_len = sizeof(mu0) / sizeof(mu0[0]);

	/* Set inner radius of the disk */
	isco = find_isco();
	printf("isco = %.15Le\n", isco);

	/* Calculate radiative efficiency */
	eta = 1.0 - specific_energy(isco);
	printf("eta = %.15Le\n", eta);

	for (int jj = 0; jj < mu_len; jj++)
	{
		cosinc = mu0[jj];
		inc = acos(cosinc); // inclination angle of the observer in rad
		printf("spin = %.15Le, inc = %.15Le, defpar = %.15Le\n", spin, inc, defpar);

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

		/* Open output file */
		sprintf(filename_o, "photons/photons4trf_a%.05Le.i%.02Le.Mdl_%.02Le.dp_%.02Le.dat", spin, cosinc, Mdl, defpar);

		foutput = fopen(filename_o, "w");

		// Start timer
		clock_t start, end, mid;
		start = clock();
		mid = clock();

		/* Assign photon position in the grid */
		for (int ii = 0; ii <= imax + 1; ii++)
		{
			if (progress_check != 0)
			{
				// Update progress for every progress_check robs
				if (ii % progress_check == 0)
				{
					// Calculate expected time usage
					iteration_time = double(clock() - mid) / double(CLOCKS_PER_SEC);
					mid = clock();
					expected_time = iteration_time * (imax - ii) / double(progress_check);
					printf("rdisk = %f; expected time left: %f minutes\n", rdisk[ii], expected_time / 60.0);
				}
			}

			// Skip first n radii
			if (ii < 50)
			{
				fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", rdisk[ii], 0.0, 0.0, 0.0, 0.0);
				continue;
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

			printf("%d gmin = %.15Le, gmax = %.15Le\n", ii, gmin, gmax);

			//*_*_*_*_*_*_*_*		CALCULATING CONSTANTLY SPACED g* GRID		*_*_*_*_*_*_*//

			gerrttol = gerrtol;

			/*---------- Branch 1 ------------*/

			xyfromrphi(rscrmin, pscrmin, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", rdiskmin, gmin, xscr, yscr, cosemmin);
			printf("%d B1:MIN %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii, rdiskmin,gmin,xscr,yscr,cosemmin);

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

					pscrcur = 0.5 * (pscrlow + pscrhigh); // Set phi_screen_current to midpoint of high and low values

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

				fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le\n", rdiskcur, gcur, xscrcur, yscrcur, cosem, gminus, xscrminus, yscrminus, gplus, xscrplus, yscrplus);
				// printf("%d B1:%d %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii,j+1,rdiskcur,gcur,xscrcur,yscrcur,cosem,gminus,xscrminus,yscrminus,gplus,xscrplus,yscrplus);
			}

			xyfromrphi(rscrmax, pscrmax, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", rdiskmax, gmax, xscr, yscr, cosemmax);
			printf("%d B1:MAX %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii,rdiskmax,gmax,xscr,yscr,cosemmax);

			/*---------- Branch 2 ------------*/

			/* Note that for Branch 2 some things are done in reverse as the values are negative */

			xyfromrphi(rscrmin, pscrmin, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", rdiskmin, gmin, xscr, yscr, cosemmin);
			printf("%d B2:MIN %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii,rdiskmin,gmin,xscr,yscr,cosemmin);

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

				fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le %.15Le\n", rdiskcur, gcur, xscrcur, yscrcur, cosem, gminus, xscrminus, yscrminus, gplus, xscrplus, yscrplus);
				// printf("%d B2:%d %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii,j+1,rdiskcur,gcur,xscrcur,yscrcur,cosem,gminus,xscrminus,yscrminus,gplus,xscrplus,yscrplus);
			}

			xyfromrphi(rscrmax, pscrmax, rdisk[ii]);

			fprintf(foutput, "%.15Le %.15Le %.15Le %.15Le %.15Le 0.0 0.0 0.0 0.0 0.0 0.0\n", rdiskmax, gmax, xscr, yscr, cosemmax);
			printf("%d B2:MAX %.6Le %.6Le %.6Le %.6Le %.6Le\n",ii,rdiskmax,gmax,xscr,yscr,cosemmax);
		}

		fclose(foutput);
	}

	return 0;
}
