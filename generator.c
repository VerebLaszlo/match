/**
 * @file generator.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "generator.h"
#include <stdio.h>

void init_Binary_System(binary_System *min, binary_System *max) {
	memset(min, 0, sizeof(binary_System));
	memset(max, 0, sizeof(binary_System));
	min->incl = DBL_MIN;
	short i;
	for (i = 0; i < 2; i++) {
		max->bh[i].m = DBL_MAX;
		max->bh[i].chi_Amp = 1.;
		max->bh[i].chi[0] = 1.;
		max->bh[i].chi[1] = 1.;
		max->bh[i].chi[2] = 1.;
		max->bh[i].theta = max->bh[i].kappa = M_PI;
		max->bh[i].varphi = max->bh[i].phi = max->bh[i].psi = 2. * M_PI;
		min->bh[i].ctheta = -1.;
		max->bh[i].ctheta = 1.;
	}
	max->M = DBL_MAX;
	max->chirpM = DBL_MAX;
	max->eta = 0.25;
	max->coaPhase = 2. * M_PI;
	max->incl = 2. * M_PI; ///< \bug nem biztos
}

void print_Binary_System(binary_System *sys, FILE *stream) {
	short i;
	fprintf(stream, "System:\n");
	fprintf(stream, "M_C="PREC"M="PREC"eta="PREC"phi_C="PREC"i="PREC"d_L="PREC"\n",
			sys->chirpM, sys->M, sys->eta, sys->coaPhase, sys->incl, sys->dist);
	for (i = 0; i < 2; i++) {
		fprintf(stream, "Black Hole %d:\n", i + 1);
		fprintf(
				stream,
				"m="PREC"th="PREC"vphi="PREC"phi="PREC"kappa="PREC"psi="PREC"ctheta="PREC"chi_Amp="PREC"chi_x="PREC"chi_y="PREC"chi_z="PREC"\n",
				sys->bh[i].m, sys->bh[i].theta, sys->bh[i].varphi,
				sys->bh[i].phi, sys->bh[i].kappa, sys->bh[i].psi,
				sys->bh[i].ctheta, sys->bh[i].chi_Amp, sys->bh[i].chi[0],
				sys->bh[i].chi[1], sys->bh[i].chi[2]);
	}
	fprintf(
			stream,
			"dec="PREC"pol="PREC"right_ascension="PREC"gmst="PREC"F_+="PREC"F_x="PREC"\n",
			sys->F.dec, sys->F.pol, sys->F.alpha, sys->F.gmst,
			sys->F.F[0], sys->F.F[1]);
	fprintf(stream, "===========================\n");fflush(stream);
}

void convert_Spins(binary_System *sys, conversion_Mode_Spins mode) {
	short i;
	if (mode == FROM_CTHETA_VPHI || mode == FROM_CTHETA_PHI) {
		sys->bh[0].theta = acos(sys->bh[0].ctheta);
		sys->bh[1].theta = acos(sys->bh[1].ctheta);
		mode--;
	}
	switch (mode) {
		case FROM_XYZ:
			for (i = 0; i < 2; i++) {
				sys->bh[i].chi_Amp = sqrt(SQR(sys->bh[i].chi[0])
						+ SQR(sys->bh[i].chi[1]) + SQR(sys->bh[i].chi[2]));
				sys->bh[i].ctheta = sys->bh[i].chi[2] / sys->bh[i].chi_Amp;
				sys->bh[i].theta = acos(sys->bh[i].ctheta);
				sys->bh[i].varphi = asin(sys->bh[i].chi[1] / sys->bh[i].chi_Amp
						/ sin(sys->bh[i].theta));
				sys->bh[i].phi = sys->bh[i].varphi + M_PI / 2.; ///< kihasználva, hogy \f$\hat{L_N} x-z síkban van, álltalánosítani kell\f$
				sys->bh[i].kappa = acos(sin(sys->bh[i].theta) * sin(
						sys->bh[i].phi) * sin(sys->incl)
						+ cos(sys->bh[i].theta) * cos(sys->incl));
				sys->bh[i].psi = atan(tan(sys->bh[i].phi) * cos(sys->incl)
						- sin(sys->incl) / (cos(sys->bh[i].phi) * tan(
								sys->bh[i].theta)));
			}
			break;
		case FROM_THETA_VPHI:
			for (i = 0; i < 2; i++) {
				sys->bh[i].ctheta = cos(sys->bh[i].theta);
				sys->bh[i].chi[0] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* cos(sys->bh[i].varphi);
				sys->bh[i].chi[1] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* sin(sys->bh[i].varphi);
				sys->bh[i].chi[2] = sys->bh[0].chi_Amp * cos(sys->bh[i].theta);
				sys->bh[i].phi = sys->bh[i].varphi + M_PI / 2.; ///< kihasználva, hogy \f$\hat{L_N} x-z síkban van, álltalánosítani kell\f$
				sys->bh[i].kappa = acos(sin(sys->bh[i].theta) * sin(
						sys->bh[i].phi) * sin(sys->incl)
						+ cos(sys->bh[i].theta) * cos(sys->incl));
				sys->bh[i].psi = atan(tan(sys->bh[i].phi) * cos(sys->incl)
						- sin(sys->incl) / (cos(sys->bh[i].phi) * tan(
								sys->bh[i].theta)));
			}
			break;
		case FROM_THETA_PHI:
			for (i = 0; i < 2; i++) {
				sys->bh[i].ctheta = cos(sys->bh[i].theta);
				sys->bh[i].kappa = acos(sin(sys->bh[i].theta) * sin(
						sys->bh[i].phi) * sin(sys->incl)
						+ cos(sys->bh[i].theta) * cos(sys->incl));
				sys->bh[i].psi = atan(tan(sys->bh[i].phi) * cos(sys->incl)
						- sin(sys->incl) / (cos(sys->bh[i].phi) * tan(
								sys->bh[i].theta)));
				sys->bh[i].varphi = sys->bh[i].phi - M_PI / 2.; ///< kihasználva, hogy \f$\hat{L_N} x-z síkban van, álltalánosítani kell\f$
				sys->bh[i].chi[0] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* cos(sys->bh[i].varphi);
				sys->bh[i].chi[1] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* sin(sys->bh[i].varphi);
				sys->bh[i].chi[2] = sys->bh[0].chi_Amp * cos(sys->bh[i].theta);
			}
			break;
		case FROM_KAPPA_PSI:
			for (i = 0; i < 2; i++) {
				sys->bh[i].ctheta = cos(sys->bh[i].kappa) * cos(sys->incl)
						- sin(sys->bh[i].kappa) * sin(sys->bh[i].phi) * sin(
								sys->incl);
				sys->bh[i].theta = acos(sys->bh[i].ctheta);
				sys->bh[i].phi = atan(tan(sys->bh[i].psi) * cos(sys->incl)
						+ sin(sys->incl) / (cos(sys->bh[i].psi) * tan(
								sys->bh[i].kappa)));
				sys->bh[i].varphi = sys->bh[i].phi - M_PI / 2.; ///< kihasználva, hogy \f$\hat{L_N} x-z síkban van, álltalánosítani kell\f$
				sys->bh[i].chi[0] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* cos(sys->bh[i].varphi);
				sys->bh[i].chi[1] = sys->bh[0].chi_Amp * sin(sys->bh[i].theta)
						* sin(sys->bh[i].varphi);
				sys->bh[i].chi[2] = sys->bh[0].chi_Amp * cos(sys->bh[i].theta);
			}
			break;
		default:
			fprintf(stderr, "Invalid spin conversion.");
			fflush(stderr);
			exit(-1);
			break;
	}
}

void convert_Masses(binary_System *sys, conversion_Mode_Masses mode) {
	switch (mode) {
		case FROM_ETAM:
			sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			sys->chirpM = pow(sys->eta, 3. / 5.) * sys->M;
			break;
		case FROM_M1M2:
			sys->M = sys->bh[0].m + sys->bh[1].m;
			sys->eta = sys->bh[0].m * sys->bh[1].m / SQR(sys->M);
			sys->chirpM = pow(sys->eta, 3. / 5.) * sys->M;
			break;
		case FROM_ETACHIRP:
			sys->M = sys->chirpM / pow(sys->eta, 3. / 5.);
			sys->bh[0].m = (1. + sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			sys->bh[1].m = (1. - sqrt(1. - 4. * sys->eta)) * sys->M / 2.;
			break;
		default:
			fprintf(stderr, "Invalid mass conversion.");
			fflush(stderr);
			exit(-1);
			break;
	}
}

// generator functions

void gen_Mass(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Masses mode) {
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

void gen_Chi(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Spin mode) {
	short i;
	switch (mode) {
		case XYZ:
			fprintf(stderr, "Not implemented yet.");
			fflush(stderr);
			exit(-1);
			break;
		case THETA_VPHI:
			for (i = 0; i < 3; i++) {
				sys->bh[i].chi_Amp
						= RANDNK(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
				sys->bh[i].theta = acos(RANDNK(cos(min->bh[i].theta), cos(max->bh[i].theta)));
				sys->bh[i].varphi
						= RANDNK(min->bh[i].varphi, max->bh[i].varphi);
			}
			convert_Spins(sys, FROM_THETA_VPHI);
			break;
		case THETA_PHI:
			for (i = 0; i < 3; i++) {
				sys->bh[i].chi_Amp
						= RANDNK(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
				sys->bh[i].theta = acos(RANDNK(cos(min->bh[i].theta), cos(max->bh[i].theta)));
				sys->bh[i].phi = RANDNK(min->bh[i].phi, max->bh[i].phi);
			}
			convert_Spins(sys, FROM_THETA_PHI);
			break;
		case KAPPA_PSI:
			for (i = 0; i < 3; i++) {
				sys->bh[i].chi_Amp
						= RANDNK(min->bh[i].chi_Amp, max->bh[i].chi_Amp);
				sys->bh[i].kappa = acos(RANDNK(cos(min->bh[i].kappa), cos(max->bh[i].kappa)));
				sys->bh[i].psi = RANDNK(min->bh[i].psi, max->bh[i].psi);
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

void gen_Sys(binary_System *sys, binary_System *min, binary_System *max) {
	sys->dist = RANDNK(min->dist, max->dist);
	sys->coaPhase = RANDNK(min->coaPhase, max->coaPhase);
	sys->incl = RANDNK(min->incl, max->incl);
	sys->F.dec = RANDNK(min->F.dec, max->F.dec);
	sys->F.pol = RANDNK(min->F.pol, max->F.pol);
	sys->F.alpha = RANDNK(min->F.alpha, max->F.alpha);
	calc_Response_For_Detector(LH, sys);
}

void gen_Parameters(binary_System *sys, binary_System *min, binary_System *max,
		gen_Mode_Masses mass, gen_Mode_Spin spin) {
	gen_Sys(sys, min, max);
	gen_Mass(sys, min, max, mass);
	gen_Chi(sys, min, max, spin);
}

void calc_Response_For_Detector(detector det, binary_System *sys) {
	double rm[3][3];
	detector_table dettable;
	double phi = sys->F.alpha - sys->F.gmst;
	dettable = GetDetectorTable(det);
	calc_Response_Matrix(dettable.nx, dettable.ny, rm);
	calc_Response(rm, sys->F.dec, phi, sys->F.pol, &(sys->F.F[0]), &(sys->F.F[1]));
}
