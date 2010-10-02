/**
 * @file generator.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "generator.h"
#include <stdio.h>

void convert_Angles_Components(binary_System *sys, conversion_Mode mode) {
	switch (mode) {
		case ANGLE_TO_COMP:
			// first spin
			sys->bh[0].chi[0] = sys->bh[0].chi_Amp * sqrt(1.
					- SQR(sys->bh[0].cth)) * cos(sys->bh[0].phi);
			sys->bh[0].chi[1] = sys->bh[0].chi_Amp * sqrt(1.
					- SQR(sys->bh[0].cth)) * sin(sys->bh[0].phi);
			sys->bh[0].chi[2] = sys->bh[0].chi_Amp * sys->bh[0].cth;
			// second spin
			sys->bh[1].chi[0] = sys->bh[1].chi_Amp * sqrt(1.
					- SQR(sys->bh[1].cth)) * cos(sys->bh[1].phi);
			sys->bh[1].chi[1] = sys->bh[1].chi_Amp * sqrt(1.
					- SQR(sys->bh[1].cth)) * sin(sys->bh[1].phi);
			sys->bh[1].chi[2] = sys->bh[1].chi_Amp * sys->bh[1].cth;
			break;
		case COMP_TO_ANGLE:
			// first spin
			sys->bh[0].chi_Amp = sqrt(SQR(sys->bh[0].chi[0])
					+ SQR(sys->bh[0].chi[1]) + SQR(sys->bh[0].chi[2]));
			sys->bh[0].cth = sys->bh[0].chi[2] / sys->bh[0].chi_Amp;
			sys->bh[0].phi = sys->bh[0].chi[1] / sys->bh[0].chi_Amp / sqrt(1.
					- SQR(sys->bh[0].cth));
			// second spin
			sys->bh[1].chi_Amp = sqrt(SQR(sys->bh[1].chi[0])
					+ SQR(sys->bh[1].chi[1]) + SQR(sys->bh[1].chi[2]));
			sys->bh[1].cth = sys->bh[1].chi[2] / sys->bh[1].chi_Amp;
			sys->bh[1].phi = sys->bh[1].chi[1] / sys->bh[1].chi_Amp / sqrt(1.
					- SQR(sys->bh[1].cth));
			break;
	}
}

void convert_etaM_m1m2(binary_System *sys, conversion_Mode mode) {
	switch (mode) {
		case ETAM_TO_M1M2:
			sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			break;
		case M1M2_TO_ETAM:
			sys->M = sys->bh[0].m + sys->bh[1].m;
			sys->eta = sys->bh[0].m * sys->bh[1].m / SQR(sys->M);
			break;
	}
}

// generator functions

binary_System lower = {
		{
				{ 0., -1., 0., 0., { 0., 0., 0. } },	// BH1
				{ 0., -1., 0., 0., { 0., 0., 0. } }		// BH2
		},
		0., 0., DBL_MIN, 1., { -M_PI_2, 0., 0., {0., 0.} } };
binary_System upper = {
		{
				{ 40., 1., 2. * M_PI, 1., { 1., 1., 1. } },
				{ 40., 1., 2. * M_PI, 1., { 1., 1., 1. } }
		},
		40., 0.25, 2. * M_PI, 10., { M_PI_2, M_PI, 2. * M_PI, {0., 0.} } };

void check_Borders(binary_System *min, binary_System *max) {
	short i;
	for (i = 0; i < 2; i++) {
		if (min->bh[i].m < lower.bh[i].m) {
			min->bh[i].m = lower.bh[i].m;
		}
		if (min->bh[i].cth < lower.bh[i].cth) {
			min->bh[i].cth = lower.bh[i].cth;
		}
		if (min->bh[i].phi < lower.bh[i].phi) {
			min->bh[i].phi = lower.bh[i].phi;
		}
		if (min->bh[i].chi_Amp < lower.bh[i].chi_Amp) {
			min->bh[i].chi_Amp = lower.bh[i].chi_Amp;
		}
		if (min->bh[i].chi[0] < lower.bh[i].chi[0]) {
			min->bh[i].chi[0] = lower.bh[i].chi[0];
		}
		if (min->bh[i].chi[1] < lower.bh[i].chi[1]) {
			min->bh[i].chi[1] = lower.bh[i].chi[1];
		}
		if (min->bh[i].chi[2] < lower.bh[i].chi[2]) {
			min->bh[i].chi[2] = lower.bh[i].chi[2];
		}
		if (max->bh[i].m > upper.bh[i].m) {
			max->bh[i].m = upper.bh[i].m;
		}
		if (max->bh[i].cth > upper.bh[i].cth) {
			max->bh[i].cth = upper.bh[i].cth;
		}
		if (max->bh[i].phi > upper.bh[i].phi) {
			max->bh[i].phi = upper.bh[i].phi;
		}
		if (max->bh[i].chi_Amp > upper.bh[i].chi_Amp) {
			max->bh[i].chi_Amp = upper.bh[i].chi_Amp;
		}
		if (max->bh[i].chi[0] > upper.bh[i].chi[0]) {
			max->bh[i].chi[0] = upper.bh[i].chi[0];
		}
		if (max->bh[i].chi[1] > upper.bh[i].chi[1]) {
			max->bh[i].chi[1] = upper.bh[i].chi[1];
		}
		if (max->bh[i].chi[2] > upper.bh[i].chi[2]) {
			max->bh[i].chi[2] = upper.bh[i].chi[2];
		}
	}
	if (min->M < lower.M) {
		min->M = lower.M;
	}
	if (max->M > upper.M) {
		max->M = upper.M;
	}
	if (min->eta < lower.eta) {
		min->eta = lower.eta;
	}
	if (max->eta > upper.eta) {
		max->eta = upper.eta;
	}
	if (min->incl < lower.incl) {
		min->incl = lower.incl;
	}
	if (max->incl > upper.incl) {
		max->incl = upper.incl;
	}
	if (min->dist < lower.dist) {
		min->dist = lower.dist;
	}
	if (max->dist > upper.dist) {
		max->dist = upper.dist;
	}
	if (min->F.dec < lower.F.dec) {
		min->F.dec = lower.F.dec;
	}
	if (min->F.pol < lower.F.pol) {
		min->F.pol = lower.F.pol;
	}
	if (min->F.phi < lower.F.phi) {
		min->F.phi = lower.F.phi;
	}
	if (max->F.dec > upper.F.dec) {
		max->F.dec = upper.F.dec;
	}
	if (max->F.pol > upper.F.pol) {
		max->F.pol = upper.F.pol;
	}
	if (max->F.phi > upper.F.phi) {
		max->F.phi = upper.F.phi;
	}
}

void gen_Mass(binary_System *sys, binary_System *min, binary_System *max,
		conversion_Mode mode) {
	switch (mode) {
		case ETAM:
			do {
				sys->M = RANDNK(min->M, max->M);
				sys->eta = RANDNK(min->eta, max->eta);
				sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
				sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			} while (min->bh[0].m > sys->bh[0].m || sys->bh[0].m > max->bh[0].m
					|| min->bh[1].m > sys->bh[1].m || sys->bh[1].m
					> max->bh[1].m);
			break;
		case M1M2:
			do {
				sys->bh[0].m = RANDNK(min->bh[0].m, max->bh[0].m);
				sys->bh[1].m = RANDNK(min->bh[1].m, max->bh[1].m);
				sys->M = sys->bh[0].m + sys->bh[1].m;
				sys->eta = sys->bh[0].m * sys->bh[1].m / SQR(sys->M);
			} while (min->M > sys->M || sys->M > max->M || min->eta > sys->eta
					|| sys->eta > max->eta);
			break;
	}
}

void gen_Chi(binary_System *sys, binary_System *min, binary_System *max) {
	//first spin
	sys->bh[0].chi_Amp = RANDNK(min->bh[0].chi_Amp, max->bh[0].chi_Amp);
	sys->bh[0].cth = RANDNK(min->bh[0].cth, max->bh[0].cth);
	sys->bh[0].phi = RANDNK(min->bh[0].phi, max->bh[0].phi);
	sys->bh[0].chi[0] = sys->bh[0].chi_Amp * sqrt(1. - SQR(sys->bh[0].cth))
			* cos(sys->bh[0].phi);
	sys->bh[0].chi[1] = sys->bh[0].chi_Amp * sqrt(1. - SQR(sys->bh[0].cth))
			* sin(sys->bh[0].phi);
	sys->bh[0].chi[2] = sys->bh[0].chi_Amp * sys->bh[0].cth;
	// second spin
	sys->bh[1].chi_Amp = RANDNK(min->bh[1].chi_Amp, max->bh[1].chi_Amp);
	sys->bh[1].cth = RANDNK(min->bh[1].cth, max->bh[1].cth);
	sys->bh[1].phi = RANDNK(min->bh[1].phi, max->bh[1].phi);
	sys->bh[1].chi[0] = sys->bh[1].chi_Amp * sqrt(1. - SQR(sys->bh[1].cth))
			* cos(sys->bh[1].phi);
	sys->bh[1].chi[1] = sys->bh[1].chi_Amp * sqrt(1. - SQR(sys->bh[1].cth))
			* sin(sys->bh[1].phi);
	sys->bh[1].chi[2] = sys->bh[1].chi_Amp * sys->bh[1].cth;
}

void gen_Sys(binary_System *sys, binary_System *min, binary_System *max) {
	sys->dist = RANDNK(min->dist, max->dist);
	sys->incl = RANDNK(min->incl + DBL_MIN, max->incl);
	sys->F.dec = RANDNK(min->F.dec, max->F.dec);
	sys->F.pol = RANDNK(min->F.dec, max->F.pol);
	sys->F.phi = RANDNK(min->F.dec, max->F.phi);
	sys->F.F[0] = sys->F.F[1] = 1. / 0.;
}

void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
		conversion_Mode mode) {
	gen_Mass(sys, min, max, mode);
	gen_Chi(sys, min, max);
	gen_Sys(sys, min, max);
}
