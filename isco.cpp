/* Code to calculate isco.dat used in python script */
// Modify deformation parameters, name of file, and upper and lower bounds as needed

#include "def.h"

int main(int argc, char *argv[])
{
	// Set deformation parameters
	defpar = 0.005;

	/* ----- Calculate ISCO for spin vs. deformation parameter space ----- */

	// Spin values
	double a[30] = {-0.998, -0.8, -0.6, -0.45, -0.3, -0.15, 0.0, 0.15, 0.3, 0.45, 0.5804822, 0.6381477, 0.6800805, 0.7142378, 0.7435757, 0.7695684, 0.7930734, 0.8146393, 0.8346417, 0.8533508, 0.8709681, 0.8876485, 0.9035142, 0.9186634, 0.9331764, 0.9471198, 0.9605495, 0.9735132, 0.9860516, 0.9982};
	// double a[1] = {0.1};
	double defparmin, defparmax, ddefpar;
	double lowbound, highbound;
	FILE *iscoout;

	int len_a = sizeof(a) / sizeof(a[0]);

	iscoout = fopen("isco.dat", "w");

	for (int i = 0; i < len_a; i++)
	{
		spin = a[i];

		//highbound = pow(1.0+sqrt(1.0-spin*spin),4.0)/(spin*spin); //Upper bound if needed
		lowbound = -2.0/3.0*1.0/2.0*pow(1.0+sqrt(1.0-spin*spin),4.0); //Lower bound if needed

		//Set upper and lower bounds if cutting off at some value
		defparmin = max(-5.0, lowbound);
		//defparmin = -5.0;
		//defparmax = min(5.0l, highbound);
		defparmax = 5.0;

		// Find isco for Kerr
		defpar = 0.0;
		isco = find_isco();
		fprintf(iscoout, "%.8Lf %.8Lf %.8Lf\n", spin, defpar, isco);

		// Find isco at max defpar value
		defpar = defparmax;
		isco = find_isco();
		fprintf(iscoout, "%.8Lf %.8Lf %.8Lf\n", spin, defpar, isco);

		// Find isco at min defpar value
		defpar = defparmin;
		isco = find_isco();
		fprintf(iscoout, "%.8Lf %.8Lf %.8Lf\n", spin, defpar, isco);

		ddefpar = (defparmax - defparmin) / 28.0;

		// Find isco for remaining values with equal spacing in defpar
		for (int j = 1, k = 1; j < 28; k++, j++)
		{
			defpar = ddefpar * k;
			isco = find_isco();
			fprintf(iscoout, "%.8Lf %.8Lf %.8Lf\n", spin, defpar, isco);

			defpar = -ddefpar * k;
			if (defpar > defparmin)
			{
				isco = find_isco();
				fprintf(iscoout, "%.8Lf %.8Lf %.8Lf\n", spin, defpar, isco);
				j++;
			}
		}
	}

	fclose(iscoout);
	return 0;
}
