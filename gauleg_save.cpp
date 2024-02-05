#include "def.h"

int main()
{
	double rdisk_f = 1000, rdisk[imax+2];
	spin = 0.1;
	isco = find_isco();
	printf("isco = %.15Le\n", isco);

	gauleg(isco, rdisk_f, rdisk);

	// Save radial grid to file
	FILE *f = fopen("gauleg.dat", "w");
	for (int i = 0; i < imax+2; i++)
		fprintf(f, "%e\n", rdisk[i]);
	fclose(f);
}