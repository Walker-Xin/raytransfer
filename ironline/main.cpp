// TODO: change confusing x->r, y->theta, z->phi
// TODO: change ODE to new form
// NOTE: changing from double to double causes some changes in the output

#include "def.h"

int main()
{
	double spin, spin2, epsilon_r, epsilon_t;
	double D, D2;
	double E_line, N_0, N_tot, N_tot1, N_tot2, alpha;
	double iobs, iobs_deg, dobs;
	double robs, pobs;
	double robs_i, robs_f, rstep, rstep2, pstep;
	double xobs, yobs, xobs2, yobs2;
	double xin, xout;
	double t0, x0, y0, phi0, r0, th0;
	double kt0, kx0, ky0, kphi0, kr0, kth0;
	double r02, s0, s02;
	double omega;
	double fact1, fact2, fact3, fact4, fact5, B, C;
	double t, x, y, phi;
	double kt, kx, ky, kphi;
	double tau, xau, yau, phiau, zau, kyau, ktau;
	double rmid, thmid;
	double kyem;
	double xmap, ymap;
	double const0, const1, b;
	double h, hstart, hnext;
	double temp_1, temp_2;
	double atol, rtol;
	double check;
	double xdist;
	double dd, ss, ssss, horizon;
	double deltax, deltay, dxdy;
	double pp, qq;
	double fr;
	double Upsilon;

	double isco;
	double g[4][4], christ[4][4][4], compare[4][4][4];
	double v[4], p[4];
	double u[4];
	// double diffs[5], vars[5], vars_temp[5], vars_4th[5], vars_5th[5], k1[5], k2[5], k3[5], k4[5], k5[5], k6[5];
	double RK1[4], RK2[4], RK3[4], RK4[4];
	double verr[4], vtol[4], err[4];
	double xem[4];
	double gfactor;
	double limbdark;

	double E_obs[imax];
	double N_obs[imax];
	double N_obs1[imax], N_obs2[imax];
	double fphi[imax], fphi0[imax];
	double fphi1[imax], fphi01[imax];
	double fphi2[imax], fphi02[imax];

	double errmin, errmax;
	double cross_tol;

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

	spin = 0.5;
	double defpar = 0;
	{
		iobs_deg = 45;
		{
			epsilon_t = 0;
			{
				epsilon_r = 0; /* deformation parameter */

				spin2 = spin * spin;

				iobs = Pi / 180 * iobs_deg; /* inclination angle of the observer in rad */

				D = 10.0; /* distance Earth-binary system in kpc */
				D2 = D * D;

				/* ----- Set model for the spectral line ----- */

				E_line = 6.4; /* energy rest of the line in keV */

				N_0 = 1.0;	/* normalization */
				alpha = -3; /* radial power law index */

				/* ----- Set inner and outer radius of the disk ----- */

				isco = find_isco();

				xin = isco;			 /* inner radius of the accretion disk; set isco */
				xin = 4.2330000000134405; // set maunally
				xout = isco + 250.0; /* outer radius of the accretion disk */

				/* ----- Set computational parameters ----- */

				dobs = 1000000; /* distance of the observer */

				robs_i = 1;
				robs_f = 150;

				rstep = 1.01;
				rstep2 = (rstep - 1) / rstep;
				pstep = 2 * Pi / 400;

				atol = 1.0e-4;
				rtol = 1.0e-4;

				hstart = 100;

				E_obs[0] = 0.1; /* minimum photon energy detected by the observer; in keV */

				for (i = 1; i <= imax - 1; i++)
					E_obs[i] = E_obs[i - 1] + 0.1;

				for (i = 0; i <= imax - 1; i++)
					N_obs[i] = 0;

				for (i = 0; i <= imax - 1; i++)
					N_obs1[i] = 0;

				for (i = 0; i <= imax - 1; i++)
					N_obs2[i] = 0;

				sprintf(filename_o, "iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat", spin, epsilon_r, epsilon_t, iobs_deg);

				/* ----- assign photon position in the grid ----- */

				for (robs = robs_i; robs < robs_f; robs = robs * rstep)
				{

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

						/* ----- compute photon initial conditions ----- */

						r02 = xobs * xobs + yobs * yobs + dobs * dobs;
						r0 = sqrt(r02);

						fact1 = dobs * sin(iobs) - yobs * cos(iobs);
						fact2 = fact1 * fact1;
						fact3 = xobs * xobs + fact2;
						fact4 = sqrt(fact3);

						t0 = 0.0;
						x0 = r0;
						y0 = acos((yobs * sin(iobs) + dobs * cos(iobs)) / r0);
						phi0 = atan(xobs / fact1);

						s0 = sin(y0);
						s02 = s0 * s0;

						kx0 = -dobs / r0;
						ky0 = (cos(iobs) - dobs * (yobs * sin(iobs) + dobs * cos(iobs)) / r02) / fact4;
						kphi0 = xobs * sin(iobs) / fact3;
						kt0 = sqrt(kx0 * kx0 + r02 * ky0 * ky0 + r02 * s02 * kphi0 * kphi0);

						/* ----- solve geodesic equations
						 fourth-order runge-kutta-nystrom method
						 see E. Lund et al., 2009 JINST 4 P04001 ----- */

						t = t0;
						x = x0;
						y = y0;
						phi = phi0;

						kt = kt0;
						kx = kx0;
						ky = ky0;
						kphi = kphi0;

						const0 = kt0;
						const1 = r02 * s02 * kphi0 / kt0;

						stop_integration = 0;

						hnext = hstart;

						do
						{

							v[0] = t;
							v[1] = x;
							v[2] = y;
							v[3] = phi;

							p[0] = kt;
							p[1] = kx;
							p[2] = ky;
							p[3] = kphi;

							do
							{

								h = hnext;

								/* ----- compute RK1 ----- */

								temp_1 = x;
								temp_2 = y;

								// christoffel_alt(spin, spin2, epsilon_r, epsilon_t, temp_1, temp_2, christ);
								christoffel(spin, defpar, temp_1, temp_2, christ);

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
												RK1[n1] = -christ[n1][0][0] * u[0] * u[0];
											}
											else
											{
												RK1[n1] -= christ[n1][n2][n3] * u[n2] * u[n3];
											}
										}
									}
								}

								/* ----- compute RK2 ----- */

								temp_1 = x + h * kx / 2 + h * h * RK1[1] / 16;
								temp_2 = y + h * ky / 2 + h * h * RK1[2] / 16;

								// christoffel_alt(spin, spin2, epsilon_r, epsilon_t, temp_1, temp_2, christ);
								christoffel(spin, defpar, temp_1, temp_2, christ);

								for (i = 0; i <= 3; i++)
									u[i] = p[i] + h * RK1[i] / 4;

								for (n1 = 0; n1 <= 3; n1++)
								{
									for (n2 = 0; n2 <= 3; n2++)
									{
										for (n3 = 0; n3 <= 3; n3++)
										{
											if (n2 == 0 && n3 == 0)
											{
												RK2[n1] = -christ[n1][0][0] * u[0] * u[0];
											}
											else
											{
												RK2[n1] -= christ[n1][n2][n3] * u[n2] * u[n3];
											}
										}
									}
								}

								/* ----- compute RK3 ----- */

								for (i = 0; i <= 3; i++)
									u[i] = p[i] + h * RK2[i] / 4;

								for (n1 = 0; n1 <= 3; n1++)
								{
									for (n2 = 0; n2 <= 3; n2++)
									{
										for (n3 = 0; n3 <= 3; n3++)
										{
											if (n2 == 0 && n3 == 0)
											{
												RK3[n1] = -christ[n1][0][0] * u[0] * u[0];
											}
											else
											{
												RK3[n1] -= christ[n1][n2][n3] * u[n2] * u[n3];
											}
										}
									}
								}

								/* ----- compute RK4 ----- */

								temp_1 = x + h * kx + h * h * RK3[1] / 4;
								temp_2 = y + h * ky + h * h * RK3[2] / 4;

								// christoffel_alt(spin, spin2, epsilon_r, epsilon_t, temp_1, temp_2, christ);
								christoffel(spin, defpar, temp_1, temp_2, christ);

								for (i = 0; i <= 3; i++)
									u[i] = p[i] + h * RK3[i] / 2;

								for (n1 = 0; n1 <= 3; n1++)
								{
									for (n2 = 0; n2 <= 3; n2++)
									{
										for (n3 = 0; n3 <= 3; n3++)
										{
											if (n2 == 0 && n3 == 0)
											{
												RK4[n1] = -christ[n1][0][0] * u[0] * u[0];
											}
											else
											{
												RK4[n1] -= christ[n1][n2][n3] * u[n2] * u[n3];
											}
										}
									}
								}

								/* ----- local error ----- */

								for (i = 0; i <= 3; i++)
								{
									verr[i] = 0.5 * h * h * (RK1[i] - RK2[i] - RK3[i] + RK4[i]);
									verr[i] *= verr[i];
									vtol[i] = atol + fabs(v[i]) * rtol;
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

							xau = x;
							yau = y;
							phiau = phi;

							kyau = ky;

							t += h * kt + (RK1[0] + RK2[0] + RK3[0]) * h * h / 12;
							x += h * kx + (RK1[1] + RK2[1] + RK3[1]) * h * h / 12;
							y += h * ky + (RK1[2] + RK2[2] + RK3[2]) * h * h / 12;
							phi += h * kphi + (RK1[3] + RK2[3] + RK3[3]) * h * h / 12;

							kt += (RK1[0] + 2 * RK2[0] + 2 * RK3[0] + RK4[0]) * h / 12;
							kx += (RK1[1] + 2 * RK2[1] + 2 * RK3[1] + RK4[1]) * h / 12;
							ky += (RK1[2] + 2 * RK2[2] + 2 * RK3[2] + RK4[2]) * h / 12;
							kphi += (RK1[3] + 2 * RK2[3] + 2 * RK3[3] + RK4[3]) * h / 12;

							/*
							 fdat = fopen("ray.dat","a");
							 fprintf(fdat,"%e %e %e %e %e\n",t,x*sin(y)*cos(phi),x*sin(y)*sin(phi),x*cos(y),mdot);
							 fclose(fdat);
							 */

							if (cos(y) < 0.0)
							{

								intersection(xau, yau, phiau, x, y, phi, xem);

								if (xem[1] > xin && xem[1] < xout)
								{
									stop_integration = 1; /* the photon hits the disk */
								}
								else
								{
									stop_integration = 2; /* the photon misses the disk */
								}
							}

							dd = x * x - 2 * x + spin2;
							ss = x * x + spin2 * cos(y) * cos(y);
							ssss = ss * ss;

							horizon = dd + spin2 * sin(y) * sin(y) * epsilon_r * x / ssss;

							if (horizon < 0.001)
								stop_integration = 4; /* the photon crosses the horizon */

							if (x < 1)
								stop_integration = 5; /* the photon hits the singularity */

							if (x != x)
								stop_integration = 6; /* numerical problems! */

							if (t < 0)
								stop_integration = 7; /* numerical problems! */

							if (x > 1.05 * dobs)
								stop_integration = 8; /* the photon escapes to infinity */

						} while (stop_integration == 0);

						if (stop_integration == 1)
						{

							redshift(spin, spin2, epsilon_r, epsilon_t, xem[1], const0, const1, kyau, gfactor, limbdark);

							/* Upsilon = 1 for isotropic radiation;
							Upsilon = limbdark[0] for limb-darkened radiation */

							/*
							Upsilon = 1;

							Upsilon = limbdark;
							*/

							pp = gfactor * E_line;

							/* --- integration - part 1 --- */

							for (i = 0; i <= imax - 2; i++)
							{

								if (E_obs[i] < pp && E_obs[i + 1] > pp)
								{

									qq = gfactor * gfactor * gfactor * gfactor;
									qq = qq * pow(xem[1], alpha);

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
			}
		}
	}

	return 0;
}