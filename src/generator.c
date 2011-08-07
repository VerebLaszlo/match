/**
 * @file generator.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "generator.h"

/**
 * XX
 * @param sys
 * @param mode
 */
void convert_Spins(binary_System *sys, conversion_Mode_Spins mode) {
	assert(sys);
	double temp;
	double PHI = M_PI / 2.0;
	double theta1 = 0.0;
	double theta2 = sys->incl;
	double xyz[3];
	switch (mode) {
	case FROM_XYZ:
		for (short i = 0; i < 2; i++) {
			sys->bh[i].chi_Amp = sqrt(
				SQR(sys->bh[i].chi[0]) + SQR(sys->bh[i].chi[1]) + SQR(sys->bh[i].chi[2]));
			sys->bh[i].ctheta = sys->bh[i].chi[2] / sys->bh[i].chi_Amp;
			sys->bh[i].theta = acos(sys->bh[i].ctheta);
			if (sys->bh[i].chi[1] < 0) {
				sys->bh[i].varphi = -acos(
					sys->bh[i].chi[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].theta));
			} else {
				sys->bh[i].varphi = acos(
					sys->bh[i].chi[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].theta));
			}
			if (!isfinite(sys->bh[i].varphi)) {
				sys->bh[i].varphi = 0.0;
			}
			xyz[0] = sys->bh[i].chi[0] * cos(theta1) * cos(theta2)
				+ sys->bh[i].chi[1] * sin(theta1) * cos(theta2) - sys->bh[i].chi[2] * sin(theta2);
			xyz[1] = -sys->bh[i].chi[0] * sin(theta1) + sys->bh[i].chi[1] * cos(theta1);
			xyz[2] = sys->bh[i].chi[0] * cos(theta1) * sin(theta2)
				+ sys->bh[i].chi[1] * sin(theta1) * sin(theta2) + sys->bh[i].chi[2] * cos(theta2);
			sys->bh[i].kappa = acos(xyz[2] / sys->bh[i].chi_Amp);
			if (xyz[1] < 0) {
				sys->bh[i].psi = -acos(xyz[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].kappa));
			} else {
				sys->bh[i].psi = acos(xyz[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].kappa));
			}
			if (!isfinite(sys->bh[i].psi)) {
				sys->bh[i].psi = 0.0;
			}
		}
		break;
	case FROM_THETA_VPHI:
		for (short i = 0; i < 2; i++) {
			sys->bh[i].ctheta = cos(sys->bh[i].theta);
			sys->bh[i].chi[0] = sys->bh[i].chi_Amp * sin(sys->bh[i].theta) * cos(sys->bh[i].varphi);
			sys->bh[i].chi[1] = sys->bh[i].chi_Amp * sin(sys->bh[i].theta) * sin(sys->bh[i].varphi);
			sys->bh[i].chi[2] = sys->bh[i].chi_Amp * cos(sys->bh[i].theta);
			xyz[0] = sys->bh[i].chi[0] * cos(theta1) * cos(theta2)
				+ sys->bh[i].chi[1] * sin(theta1) * cos(theta2) - sys->bh[i].chi[2] * sin(theta2);
			xyz[1] = -sys->bh[i].chi[0] * sin(theta1) + sys->bh[i].chi[1] * cos(theta1);
			xyz[2] = sys->bh[i].chi[0] * cos(theta1) * sin(theta2)
				+ sys->bh[i].chi[1] * sin(theta1) * sin(theta2) + sys->bh[i].chi[2] * cos(theta2);
			sys->bh[i].kappa = acos(xyz[2] / sys->bh[i].chi_Amp);
			if (xyz[1] < 0) {
				sys->bh[i].psi = -acos(xyz[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].kappa));
			} else {
				sys->bh[i].psi = acos(xyz[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].kappa));
			}
			if (!isfinite(sys->bh[i].psi)) {
				sys->bh[i].psi = 0.0;
			}
		}
		break;
	case FROM_THETA_PHI:
		for (short i = 0; i < 2; i++) {
			sys->bh[i].varphi = sys->bh[i].phi - PHI;
			sys->bh[i].ctheta = cos(sys->bh[i].theta);
			temp = cos(sys->bh[i].theta) * cos(sys->incl) - //
				sin(sys->bh[i].theta) * cos(sys->bh[i].phi) * //
					sin(sys->incl);
			if (1. - temp < 1.e-10 || 1. + temp < 1.e-10)
				temp = 1.;
			if (temp < 1.e-10)
				temp = 0.;
			sys->bh[i].kappa = acos(temp);
			temp = (sin(sys->bh[i].theta) * cos(sys->bh[i].phi) * //
				cos(sys->incl) + cos(sys->bh[i].theta) * //
				sin(sys->incl)) / sin(sys->bh[i].kappa);
			if (fabs(temp) < 1.e-15) {
				temp = 0.;
			}
			sys->bh[i].psi = acos(temp);
			sys->bh[i].chi[0] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta) * cos(sys->bh[i].varphi);
			sys->bh[i].chi[1] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta) * sin(sys->bh[i].varphi);
			sys->bh[i].chi[2] = sys->bh[0].chi_Amp * cos(sys->bh[i].theta);
		}
		break;
	case FROM_KAPPA_PSI:
		if (theta1 != 0.0) {
			theta1 *= -1.0;
		}
		if (theta2 != 0.0) {
			theta2 *= -1.0;
		}
		for (short i = 0; i < 2; i++) {
			double cos_psi;
			if (sys->bh[i].psi == M_PI_2) {
				cos_psi = 0.0;
			} else {
				cos_psi = cos(sys->bh[i].psi);
			}
			double cos_kappa;
			if (sys->bh[i].kappa == M_PI_2) {
				cos_kappa = 0.0;
			} else {
				cos_kappa = cos(sys->bh[i].kappa);
			}
			xyz[0] = sin(sys->bh[i].kappa) * cos_psi;
			xyz[1] = sin(sys->bh[i].kappa) * sin(sys->bh[i].psi);
			xyz[2] = cos_kappa;
			sys->bh[i].chi[0] = sys->bh[i].chi_Amp
				* (xyz[0] * cos(theta1) * cos(theta2) + xyz[1] * sin(theta1)
					- xyz[2] * cos(theta1) * sin(theta2));
			sys->bh[i].chi[1] = sys->bh[i].chi_Amp
				* (-xyz[0] * sin(theta1) * sin(theta2) + xyz[1] * cos(theta1)
					+ xyz[2] * sin(theta1) * sin(theta2));
			sys->bh[i].chi[2] = sys->bh[i].chi_Amp * (xyz[0] * sin(theta2) + xyz[2] * cos(theta2));
			sys->bh[i].ctheta = sys->bh[i].chi[2] / sys->bh[i].chi_Amp;
			sys->bh[i].theta = acos(sys->bh[i].ctheta);
			if (sys->bh[i].chi[1] < 0) {
				sys->bh[i].varphi = -acos(
					sys->bh[i].chi[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].theta));
			} else {
				sys->bh[i].varphi = acos(
					sys->bh[i].chi[0] / sys->bh[i].chi_Amp / sin(sys->bh[i].theta));
			}
			if (!isfinite(sys->bh[i].varphi)) {
				sys->bh[i].varphi = 0.0;
			}
		}
		break;
	default:
		fprintf(stderr, "Invalid spin conversion.");
		fflush(stderr);
		exit(EXIT_FAILURE);
		break;
	}
}

/**
 * XX
 * @param sys
 * @param mode
 */
void convert_Masses(binary_System *sys, conversion_Mode_Masses mode) {
	assert(sys);
	switch (mode) {
	case FROM_ETAM:
		sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
		sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
		sys->chirpM = pow(sys->eta, 3. / 5.) * sys->M;
		sys->mu =
			sys->bh[0].m > sys->bh[1].m ? sys->bh[0].m / sys->bh[1].m : sys->bh[1].m / sys->bh[0].m;
		break;
	case FROM_M1M2:
		sys->M = sys->bh[0].m + sys->bh[1].m;
		sys->eta = sys->bh[0].m * sys->bh[1].m / SQR(sys->M);
		sys->chirpM = pow(sys->eta, 3. / 5.) * sys->M;
		sys->mu =
			sys->bh[0].m > sys->bh[1].m ? sys->bh[0].m / sys->bh[1].m : sys->bh[1].m / sys->bh[0].m;
		break;
	case FROM_ETACHIRP:
		sys->M = sys->chirpM / pow(sys->eta, 3. / 5.);
		sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
		sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
		sys->mu =
			sys->bh[0].m > sys->bh[1].m ? sys->bh[0].m / sys->bh[1].m : sys->bh[1].m / sys->bh[0].m;
		break;
	default:
		fprintf(stderr, "Invalid mass conversion.");
		fflush(stderr);
		exit(EXIT_FAILURE);
		break;
	}
}

// generator functions

/**
 * XX
 * @param sys
 * @param min
 * @param max
 * @param mode
 */
static void gen_Mass(binary_System *sys, binary_System *min, binary_System *max,
	gen_Mode_Masses mode) {
	assert(sys);
	assert(min);
	assert(max);
	switch (mode) {
	case ETAM:
		do {
			sys->M = randomBetween(min->M, max->M);
			sys->eta = randomBetween(min->eta, max->eta);
			sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
		} while (min->bh[0].m > sys->bh[0].m || sys->bh[0].m > max->bh[0].m
			|| min->bh[1].m > sys->bh[1].m || sys->bh[1].m > max->bh[1].m);
		sys->mu =
			sys->bh[0].m > sys->bh[1].m ? sys->bh[0].m / sys->bh[1].m : sys->bh[1].m / sys->bh[0].m;
		sys->chirpM = pow(sys->M, 2.0 / 5.0) * pow(sys->eta, 3.0 / 5.0);
		break;
	case M1M2:
		do {
			sys->bh[0].m = randomBetween(min->bh[0].m, max->bh[0].m);
			sys->bh[1].m = randomBetween(min->bh[1].m, max->bh[1].m);
			sys->M = sys->bh[0].m + sys->bh[1].m;
			sys->eta = sys->bh[0].m * sys->bh[1].m / SQR(sys->M);
		} while (min->M > sys->M || sys->M > max->M || min->eta > sys->eta || sys->eta > max->eta);
		sys->mu =
			sys->bh[0].m > sys->bh[1].m ? sys->bh[0].m / sys->bh[1].m : sys->bh[1].m / sys->bh[0].m;
		break;
	case ETACHRIP:
		fprintf(stderr, "Not implemented yet.");
		fflush(stderr);
		exit(-1);
		break;
	default:
		fprintf(stderr, "Invalid mass generation mode.");
		fflush(stderr);
		exit(-1);
		break;
	}
}

/**
 * XX
 * @param sys
 * @param min
 * @param max
 * @param mode
 */
static void gen_Chi(binary_System *sys, binary_System *min, binary_System *max, gen_Mode_Spin mode) {
	assert(sys);
	assert(min);
	assert(max);
	switch (mode) {
	case XYZ:
		fprintf(stderr, "Not implemented yet.");
		fflush(stderr);
		exit(-1);
		break;
	case THETA_VPHI:
		for (short i = 0; i < 3; i++) {
			sys->bh[i].chi_Amp = randomBetween(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
			sys->bh[i].theta = acos(randomBetween(cos(min->bh[i].theta), cos(max->bh[i].theta)));
			sys->bh[i].varphi = randomBetween(min->bh[i].varphi, max->bh[i].varphi);
		}
		convert_Spins(sys, FROM_THETA_VPHI);
		break;
	case THETA_PHI:
		for (short i = 0; i < 3; i++) {
			sys->bh[i].chi_Amp = randomBetween(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
			sys->bh[i].theta = acos(randomBetween(cos(min->bh[i].theta), cos(max->bh[i].theta)));
			sys->bh[i].phi = randomBetween(min->bh[i].phi, max->bh[i].phi);
		}
		convert_Spins(sys, FROM_THETA_PHI);
		break;
	case KAPPA_PSI:
		for (short i = 0; i < 3; i++) {
			sys->bh[i].chi_Amp = randomBetween(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
			sys->bh[i].kappa = acos(randomBetween(cos(min->bh[i].kappa), cos(max->bh[i].kappa)));
			sys->bh[i].psi = randomBetween(min->bh[i].psi, max->bh[i].psi);
		}
		convert_Spins(sys, FROM_KAPPA_PSI);
		break;
	default:
		fprintf(stderr, "Invalid spin generation mode.");
		fflush(stderr);
		exit(-1);
		break;
	}
}

/**
 * XX
 * @param sys
 * @param min
 * @param max
 */
static void gen_Sys(binary_System *sys, binary_System *min, binary_System *max) {
	assert(sys);
	assert(min);
	assert(max);
	sys->dist = randomBetween(min->dist, max->dist);
	sys->coaPhase = randomBetween(min->coaPhase, max->coaPhase);
	sys->incl = randomBetween(min->incl, max->incl);
	sys->F.declination = randomBetween(min->F.declination, max->F.declination);
	sys->F.polarization = randomBetween(min->F.polarization, max->F.polarization);
	sys->F.right_Ascention = randomBetween(min->F.right_Ascention, max->F.right_Ascention);
	sys->F.gmst = randomBetween(min->F.gmst, max->F.gmst);
	calc_Antenna_Pattern_For(LH, &(sys->F));
}

/**
 * XX
 * @param sys
 * @param min
 * @param max
 * @param mass
 * @param spin
 */
void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
	gen_Mode_Masses mass, gen_Mode_Spin spin) {
	assert(sys);
	assert(min);
	assert(max);
	gen_Sys(sys, min, max);
	gen_Mass(sys, min, max, mass);
	gen_Chi(sys, min, max, spin);
}

void read_Binary_Parameter_Limits(FILE*file, binary_System *sys) {
	short length = 100;
	char line[length];
	fgets(line, length, file);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].M, &sys[1].M);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].eta, &sys[1].eta);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[0].m, &sys[1].bh[0].m);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[1].m, &sys[1].bh[1].m);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[0].chi_Amp, &sys[1].bh[0].chi_Amp);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[1].chi_Amp, &sys[1].bh[1].chi_Amp);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[0].kappa, &sys[1].bh[0].kappa);
	sys[0].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[1].kappa, &sys[1].bh[1].kappa);
	sys[0].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[0].psi, &sys[1].bh[0].psi);
	sys[0].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].bh[1].psi, &sys[1].bh[1].psi);
	sys[0].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].incl, &sys[1].incl);
	sys[0].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].dist, &sys[1].dist);
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].F.polarization, &sys[1].F.polarization);
	sys[0].F.polarization *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].F.polarization *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].F.right_Ascention, &sys[1].F.right_Ascention);
	sys[0].F.right_Ascention *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	sys[1].F.right_Ascention *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].F.declination, &sys[1].F.declination);
	sys[0].F.declination *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	sys[1].F.declination *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fgets(line, length, file);
	sscanf(line, "%lg %lg %*s\n", &sys[0].F.gmst, &sys[1].F.gmst);
	sys[0].F.declination *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	sys[1].F.declination *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
}

void print_Binary_Parameter_Limits(FILE* file, binary_System*sys) {
	fprintf(file, "%10s %10.4lg %10.4lg\n", "M", sys[0].M, sys[1].M);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "eta", sys[0].eta, sys[1].eta);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "m1", sys[0].bh[0].m, sys[1].bh[0].m);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "m2", sys[0].bh[1].m, sys[1].bh[1].m);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "chi1", sys[0].bh[0].chi_Amp, sys[1].bh[0].chi_Amp);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "chi2", sys[0].bh[1].chi_Amp, sys[1].bh[1].chi_Amp);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "kappa1", sys[0].bh[0].kappa, sys[1].bh[0].kappa);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "kappa2", sys[0].bh[1].kappa, sys[1].bh[1].kappa);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "psi1", sys[0].bh[0].psi, sys[1].bh[0].psi);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "psi2", sys[0].bh[1].psi, sys[1].bh[1].psi);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "incl", sys[0].incl, sys[1].incl);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "dist", sys[0].dist, sys[1].dist);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "pol", sys[0].F.polarization, sys[1].F.polarization);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "ra", sys[0].F.right_Ascention,
		sys[1].F.right_Ascention);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "dec", sys[0].F.declination, sys[1].F.declination);
	fprintf(file, "%10s %10.4lg %10.4lg\n", "gmst", sys[0].F.gmst, sys[1].F.gmst);
}

void print_Binary_Parameter_For_Plot(FILE* file, binary_System*sys) {
	fprintf(file, "%-13s %10.4lg %10.4lg %10.4lg\n", "#masses    ", sys[0].bh[0].m, sys[0].bh[1].m,
		sys[0].bh[1].m / sys[0].bh[0].m);
	fprintf(file, "%-13s %10.4lg %10.4lg %10.4lg\n", "#spin1     ", sys[0].bh[0].chi_Amp,
		sys[0].bh[0].kappa, sys[0].bh[0].psi);
	fprintf(file, "%-13s %10.4lg %10.4lg %10.4lg\n", "#spin2     ", sys[0].bh[1].chi_Amp,
		sys[0].bh[1].kappa, sys[0].bh[1].psi);
	fprintf(file, "%-13s %10.4lg %10.4lg\n", "#incl, dist", sys[0].incl, sys[0].dist);
}
