void diffeqs(double vars[], double diffs[])
{
	double t, r, th, phi;
	double dt, dr, dth, dphi;
	double christ[4][4][4];

	t = vars[0];
	r = vars[1];
	th = vars[2];
	phi = vars[3];

	dt = vars[4];
	dr = vars[5];
	dth = vars[6];
	dphi = vars[7];

	christoffel(r, th, christ);
	// Christoffel_jiale(r, th, christ);

	/* 1st order diff eqs */
	diffs[0] = dt;	 // dt
	diffs[1] = dr;	 // dr
	diffs[2] = dth;	 // dth
	diffs[3] = dphi; // dphi

	/* 2nd order diff eqs for r and theta; c.f. Eq (28) and Eq (29) */
	diffs[4] = -christ[0][0][0] * dt * dt - 2.0 * christ[0][0][1] * dt * dr - 2.0 * christ[0][0][2] * dt * dth - 2.0 * christ[0][0][3] * dt * dphi - christ[0][1][1] * dr * dr - 2.0 * christ[0][1][2] * dr * dth - 2.0 * christ[0][1][3] * dr * dphi - christ[0][2][2] * dth * dth - 2.0 * christ[0][2][3] * dth * dphi - christ[0][3][3] * dphi * dphi; // d2t
	diffs[5] = -christ[1][0][0] * dt * dt - 2.0 * christ[1][0][1] * dt * dr - 2.0 * christ[1][0][2] * dt * dth - 2.0 * christ[1][0][3] * dt * dphi - christ[1][1][1] * dr * dr - 2.0 * christ[1][1][2] * dr * dth - 2.0 * christ[1][1][3] * dr * dphi - christ[1][2][2] * dth * dth - 2.0 * christ[1][2][3] * dth * dphi - christ[1][3][3] * dphi * dphi; // d2r
	diffs[6] = -christ[2][0][0] * dt * dt - 2.0 * christ[2][0][1] * dt * dr - 2.0 * christ[2][0][2] * dt * dth - 2.0 * christ[2][0][3] * dt * dphi - christ[2][1][1] * dr * dr - 2.0 * christ[2][1][2] * dr * dth - 2.0 * christ[2][1][3] * dr * dphi - christ[2][2][2] * dth * dth - 2.0 * christ[2][2][3] * dth * dphi - christ[2][3][3] * dphi * dphi; // d2th
	diffs[7] = -christ[3][0][0] * dt * dt - 2.0 * christ[3][0][1] * dt * dr - 2.0 * christ[3][0][2] * dt * dth - 2.0 * christ[3][0][3] * dt * dphi - christ[3][1][1] * dr * dr - 2.0 * christ[3][1][2] * dr * dth - 2.0 * christ[3][1][3] * dr * dphi - christ[3][2][2] * dth * dth - 2.0 * christ[3][2][3] * dth * dphi - christ[3][3][3] * dphi * dphi; // d2phi
}
