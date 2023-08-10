/* Calculate cosine of the emission angle modulo the redshift factor */
long double emis_angle(long double r, long double th, long double kr, long double kth)
{
	long double Zr, Zth, k, cs1, ss1, ss3, r3, angle;
	long double gupper[4][4];

	cs1 = cos(th);
	ss1 = sin(th);
	ss3 = ss1 * ss1 * ss1;
	r3 = r * r * r;

	k = 3.0 / eta * Mdl;
	Zr = 0.5 * k * sqrt(isco / (r3 * ss1)) - cs1;
	Zth = 0.5 * k * cs1 * sqrt(isco / (r * ss3)) + r * ss1;

	uppermetric(r, th, gupper);

	angle = 1.0 / sqrt(gupper[1][1] * Zr * Zr + gupper[2][2] * Zth * Zth) * (Zr * kr + Zth * kth); // c.f. Eq (16) in Abdiakamalov 2020

	if (angle < 0.0)
		angle *= -1.0;

	return angle;
}