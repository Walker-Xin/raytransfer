// TODO: change RK45 to RKN method
// NOTE: changing from double to long double causes some changes in the output

#include "def.h"

int main(int argc, char *argv[])
{
	double D, D2;
	double E_line, N_0, N_tot, alpha;
	double robs, pobs;
	double robs_i, robs_f, rstep, rstep2, pstep;
	double xobs, yobs;
	double rin, rout;
	double pp, qq;
	double fr;

	double traced[2];
	double gfactor;

	double E_obs[imax];
	double N_obs[imax];
	double N_obs1[imax], N_obs2[imax];
	double fphi[imax], fphi0[imax];
	double fphi1[imax], fphi01[imax];
	double fphi2[imax], fphi02[imax];

	double errtol, errmin, errmax;

	int errcheck, crosscheck = 0, acccheck = 0, blockcheck = 0;
	int stop_integration;
	int n1, n2, n3;
	int i, j, m;
	int n_robs, i_robs;

	char filename_i[128];
	char filename_o[128];

	double time_taken, iteration_time, expected_time;

	FILE *finput;
	FILE *foutput;
	FILE *ftab;
	FILE *fdat;
	FILE *fgoal;

	// Set default computation parameters
	spin = 0.5; // spin parameter
	spin2 = spin * spin;
	defpar = 0.5;	 // deformation parameter
	iobs_deg = 45.0; // inclination angle of the observer in degrees
	dobs = 1.0e+8;
	errtol = 1.0e-8; // error tolerance for RK45

	// Set computation parameters from user input if provided
	if (argc > 1)
	{
		spin = atof(argv[1]); // spin parameter
		defpar = atof(argv[2]);
		iobs_deg = atof(argv[3]); // inclination angle of the observer in degrees
		errtol = atof(argv[4]);	  // error tolerance for RK45
		printf("Using user input parameters. spin = %f, deformation = %f, inclination = %f, error tolerence = %e\n", spin, defpar, iobs_deg, errtol);
	}
	else
	{
		printf("Using default parameters. spin = %f, deformation = %f, inclination = %f, error tolerence = %e\n", spin, defpar, iobs_deg, errtol);
	}

	inc = Pi / 180 * iobs_deg; // inclination angle of the observer in rad
	D = 10.0;				   /* distance Earth-binary system in kpc */
	D2 = D * D;

	// Set model for the spectral line
	E_line = 6.4; // energy rest of the line in keV

	N_0 = 1.0;	// normalization
	alpha = -3; // radial power law index

	// Set inner and outer radius of the disk
	isco = find_isco(spin, defpar);
	// isco = 4.2330000000134405; // set maunally
	printf("Innermost stable circular orbit: %f\n", isco);

	rin = isco;			 /* inner radius of the accretion disk; set isco */
	rout = isco + 250.0; /* outer radius of the accretion disk */

	/* ----- Set computational parameters ----- */

	robs_i = 1.0;
	robs_f = 150.0;

	rstep = 1.01;
	rstep2 = (rstep - 1) / rstep;
	pstep = 2 * Pi / 400;

	errmin = errtol / 10.0;
	errmax = errtol * 10.0;

	E_obs[0] = 0.1; /* minimum photon energy detected by the observer; in keV */

	for (i = 1; i <= imax - 1; i++)
		E_obs[i] = E_obs[i - 1] + 0.1;

	for (i = 0; i <= imax - 1; i++)
		N_obs[i] = 0;

	for (i = 0; i <= imax - 1; i++)
		N_obs1[i] = 0;

	for (i = 0; i <= imax - 1; i++)
		N_obs2[i] = 0;

	sprintf(filename_o, "iron_a%.03f.def%.02f.i%.02f.dat", spin, defpar, iobs_deg);

	// Start timer
	clock_t start, end, mid;
	start = clock();
	mid = clock();

	// Calculate number of robs
	n_robs = 0;
	i_robs = 0;
	for (robs = robs_i; robs < robs_f; robs = robs * rstep)
	{
		n_robs++;
	}

	/* ----- assign photon position in the grid ----- */

	for (robs = robs_i; robs < robs_f; robs = robs * rstep)
	{
		i_robs++;

		// Update progress for every 10 robs
		if (i_robs % 10 == 0)
		{
			// Calculate expected time usage
			iteration_time = double(clock() - mid) / double(CLOCKS_PER_SEC);
			mid = clock();
			expected_time = iteration_time * (n_robs - i_robs) / 10.0;
			printf("robs = %f; expected time left: %f minutes\n", robs, expected_time / 60.0);
		}

		for (i = 0; i <= imax - 1; i++)
			fphi[i] = 0;
		for (i = 0; i <= imax - 1; i++)
			fphi1[i] = 0;
		for (i = 0; i <= imax - 1; i++)
			fphi2[i] = 0;

		for (pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep)
		{
			xobs = robs * cos(pobs);
			yobs = robs * sin(pobs) * cos(inc);

			stop_integration = raytrace(errmin, errmax, xobs, yobs, traced);

			if (stop_integration == 1 && traced[1] != 0)
			{

				gfactor = traced[1];
				// printf("%e %e\n", traced[0], gfactor);

				pp = gfactor * E_line;

				/* --- integration - part 1 --- */

				for (i = 0; i <= imax - 2; i++)
				{

					if (E_obs[i] < pp && E_obs[i + 1] > pp)
					{

						qq = gfactor * gfactor * gfactor * gfactor;
						qq = qq * pow(traced[2], alpha);

						fphi[i] = fphi[i] + qq;
					}
				}
			}
		}

		/* --- integrazion - part 2 --- */

		for (i = 0; i <= imax - 1; i++)
		{

			fr = robs * robs * fphi[i] * rstep2;

			N_obs[i] = N_obs[i] + fr;
		}
	}

	// End timer and print time taken in miniutes
	end = clock();
	time_taken = double(end - start) / double(CLOCKS_PER_SEC);
	printf("Time taken by program is : %f minutes\n", time_taken / 60.0);

	/* --- print spectrum --- */

	foutput = fopen(filename_o, "w");

	N_tot = 0.0;

	for (i = 0; i <= imax - 1; i++)
	{

		N_obs[i] = N_0 * N_obs[i] * pstep * cos(inc) / E_obs[i];

		N_tot = N_tot + N_obs[i];
	}

	for (i = 0; i <= imax - 1; i++)
	{

		fprintf(foutput, "%f %f\n", E_obs[i], N_obs[i] / N_tot);
	}

	fclose(foutput);

	return 0;
}