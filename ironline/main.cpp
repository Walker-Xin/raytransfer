// TODO: change RK45 to RKN method
// NOTE: changing from double to double causes some changes in the output

#include "def.h"

int main()
{
	double spin, spin2, epsilon_r, epsilon_t;
	double D, D2;
	double E_line, N_0, N_tot, alpha;
	double iobs, iobs_deg, dobs;
	double robs, pobs;
	double robs_i, robs_f, rstep, rstep2, pstep;
	double xobs, yobs;
	double rin, rout;
	double pp, qq;
	double fr;

	double isco;
	double traced[2];
	double gfactor;

	double E_obs[imax];
	double N_obs[imax];
	double N_obs1[imax], N_obs2[imax];
	double fphi[imax], fphi0[imax];
	double fphi1[imax], fphi01[imax];
	double fphi2[imax], fphi02[imax];

	double errmin, errmax;

	int errcheck, crosscheck = 0, acccheck = 0, blockcheck = 0;
	int stop_integration;
	int n1, n2, n3;
	int i, j, m;

	char filename_i[128];
	char filename_o[128];

	FILE *finput;
	FILE *foutput;
	FILE *ftab;
	FILE *fdat;
	FILE *fgoal;

	/* ----- Set free parameters ----- */

	spin = 0.9;
	defpar = 0;
	iobs_deg = 45;
	epsilon_t = 0;
	epsilon_r = 0; /* deformation parameter */

	spin2 = spin * spin;

	iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */
	inc = iobs;

	D = 10.0; /* distance Earth-binary system in kpc */
	D2 = D * D;

	/* ----- Set model for the spectral line ----- */

	E_line = 6.4; /* energy rest of the line in keV */

	N_0 = 1.0;	/* normalization */
	alpha = -3; /* radial power law index */

	/* ----- Set inner and outer radius of the disk ----- */

	isco = find_isco(spin, defpar);
	// isco = 4.2330000000134405; // set maunally

	rin = isco;			 /* inner radius of the accretion disk; set isco */
	rout = isco + 250.0; /* outer radius of the accretion disk */

	/* ----- Set computational parameters ----- */

	dobs = 1000000; /* distance of the observer */

	robs_i = 1;
	robs_f = 150;

	rstep = 1.01;
	rstep2 = (rstep - 1) / rstep;
	pstep = 2 * Pi / 100;

    errmin = 1.0e-9; // error bounds for RK45
    errmax = 1.0e-7;

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

	/* ----- assign photon position in the grid ----- */

	for (robs = robs_i; robs < robs_f; robs = robs * rstep)
	{
		printf("Now at robs = %e\n", robs);

		for (i = 0; i <= imax - 1; i++)
			fphi[i] = 0;
		for (i = 0; i <= imax - 1; i++)
			fphi1[i] = 0;
		for (i = 0; i <= imax - 1; i++)
			fphi2[i] = 0;

		for (pobs = 0; pobs < 2 * Pi - 0.5 * pstep; pobs = pobs + pstep)
		{
			xobs = robs * cos(pobs);
			yobs = robs * sin(pobs) * cos(iobs);

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

	/* --- print spectrum --- */

	foutput = fopen(filename_o, "w");

	N_tot = 0.0;

	for (i = 0; i <= imax - 1; i++)
	{

		N_obs[i] = N_0 * N_obs[i] * pstep * cos(iobs) / E_obs[i];

		N_tot = N_tot + N_obs[i];
	}

	for (i = 0; i <= imax - 1; i++)
	{

		fprintf(foutput, "%e %e\n", E_obs[i], N_obs[i] / N_tot);
	}

	fclose(foutput);

	return 0;
}