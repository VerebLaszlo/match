/**
 * @file generator.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "generator.h"

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
