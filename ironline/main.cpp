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
	char filename_tab[128];

	double time_taken, iteration_time, expected_time;
	int progress_check, first_check = 0;
	int ode_solver = 0; // 0 for RK45, 1 for RKN

	FILE *finput;
	FILE *foutput;
	FILE *ftab;
	FILE *fdat;
	FILE *fgoal;

	// Set default computation parameters
	spin = 0.5; // spin parameter
	defpar = 0.0;	 // deformation parameter
	iobs_deg = 20.0; // inclination angle of the observer in degrees
	dobs = 1.0e+6; // distance Earth-binary system in kpc
	errtol = 1.0e-6; // error tolerance
	rstep = 1.008; // step size for robs
	pstep = 2 * Pi / 800; // step size for pobs
	progress_check = 1; // check progress for every 20 robs
	ode_solver = 0; // 0 for RK45, 1 for RKN, 2 for RKN_bambi

	// Set computation parameters from user input if provided
	if (argc > 1)
	{
		spin = atof(argv[1]); // spin parameter
		defpar = atof(argv[2]);
		iobs_deg = atof(argv[3]); // inclination angle of the observer in degrees
		errtol = atof(argv[4]);	  // error tolerance for RK45
		rstep = atof(argv[5]);
		pstep = atof(argv[6]);
		progress_check = atoi(argv[7]);
		ode_solver = atoi(argv[8]);
		printf("Using user input parameters. spin=%f, deformation=%f, inclination=%f, error tolerance=%e\n, rstep=%f, pstep=%f\n", spin, defpar, iobs_deg, errtol, rstep, pstep);
	}
	else
	{
		printf("Using preset parameters. spin=%f, deformation=%f, inclination=%f, error tolerance=%e\n, rstep=%f, pstep=%f\n", spin, defpar, iobs_deg, errtol, rstep, pstep);
	}

	if (progress_check != 0)
	{
		printf("Progress check for every %d robs\n", progress_check);
	}

	if (ode_solver == 0)
	{
		printf("--------------------------------------------------\n");
		printf("Using RK45 method\n");
		printf("--------------------------------------------------\n");
	}
	else if (ode_solver == 1)
	{
		printf("--------------------------------------------------\n");
		printf("Using RKN method\n");
		printf("--------------------------------------------------\n");
	}
	else if (ode_solver == 2)
	{
		printf("--------------------------------------------------\n");
		printf("Using RKN_bambi method\n");
		printf("--------------------------------------------------\n");
	}
	else
	{
		printf("Invalid ode_solver input. Aborting.\n");
		exit(0);
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
	// isco = 4.2330000000134405; // set manually
	printf("Innermost stable circular orbit: %f\n", isco);

	rin = isco;			 /* inner radius of the accretion disk; set isco */
	rout = isco + 250.0; /* outer radius of the accretion disk */

	horizon = 1. + sqrt(1. - spin2); // horizon radius

	// Set additional computational parameters
	robs_i = 1.0;
	robs_f = 150.0;
	spin2 = spin * spin;
	rstep2 = (rstep - 1) / rstep;
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

	if (ode_solver == 0)
	{
		sprintf(filename_o, "iron_a%.03f.def%.02f.i%.02f.dat", spin, defpar, iobs_deg);
	}
	else if (ode_solver == 1)
	{
		sprintf(filename_o, "iron_a%.03f.def%.02f.i%.02f.RKN.dat", spin, defpar, iobs_deg);
	}
	else if (ode_solver == 2)
	{
		sprintf(filename_o, "iron_a%.03f.def%.02f.i%.02f.RKNb.dat", spin, defpar, iobs_deg);
	}
	else
	{
		printf("Invalid ode_solver input. Aborting.\n");
		exit(0);
	}

	// Create output observation file
	sprintf(filename_tab, "observed_a%.03f.def%.02f.i%.02f.dat", spin, defpar, iobs_deg);

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
		if (progress_check != 0)
		{
			i_robs++;

			// Update progress for every progress_check robs
			if (i_robs % progress_check == 0)
			{
				// Calculate expected time usage
				iteration_time = double(clock() - mid) / double(CLOCKS_PER_SEC);
				mid = clock();
				expected_time = iteration_time * (n_robs - i_robs) / double(progress_check);
				printf("robs = %f; expected time left: %f minutes\n", robs, expected_time / 60.0);
			}
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

			if (ode_solver == 0){
				stop_integration = raytrace(errmin, errmax, xobs, yobs, traced);
			}
			else if (ode_solver == 1){
				stop_integration = raytrace_RKN(errtol, xobs, yobs, traced);
			}
			else if (ode_solver == 2){
				stop_integration = raytrace_RKN_bambi(errtol, xobs, yobs, traced);
			}
			else{
				printf("Invalid ode_solver input. Aborting.\n");
				exit(0);
			}

			if (stop_integration == 1 && traced[1] != 0)
			{
				// Record initial photon positions
				foutput = fopen(filename_tab, "a");
				fprintf(foutput, "%f %f\n", xobs, yobs);
				fclose(foutput);

				if (first_check == 0)
				{
					printf("First check: robs = %f, pobs = %f, rmid = %f, gfactor = %f, xem[1] = %f\n", robs, pobs, traced[0], traced[1], traced[2]);
					first_check = 1;
					// exit(0);
				}

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

		/* --- integration - part 2 --- */

		for (i = 0; i <= imax - 1; i++)
		{

			fr = robs * robs * fphi[i] * rstep2;

			N_obs[i] = N_obs[i] + fr;
		}
	}

	// End timer and print time taken in minutes
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

	// Show output name
	printf("Output file: %s\n", filename_o);

	return 0;
}
