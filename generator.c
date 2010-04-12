/**
 * @file generator.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "generator.h"

binary_system min, max;

void gen_M_eta(dpc m1, dpc m2, dpc M, dpc eta) {
	do {
		*M = min.M + (max.M - min.M) * rand1();
		// \todo át kell-e állítani eta-t?
		*eta = min.eta + (max.eta - min.eta) * rand1();
		double temp = sqrt(1. - 4. * *eta);
		*m1 = (*M / 2.) * (1. + temp);
		*m2 = (*M / 2.) * (1. - temp);
	} while (min.bh1.m > *m1 || *m1 > max.bh1.m || min.bh2.m > *m2 || *m2 > max.bh2.m);
}

void gen_m1_m2(dpc m1, dpc m2, dpc M, dpc eta) {
	*m1 = min.bh1.m + (max.bh1.m - min.bh1.m) * rand1();
	*m2 = min.bh2.m + (max.bh2.m - min.bh2.m) * rand1();
	*M = *m1 + *m2;
	*eta = *m1 * *m2 / (*M * *M);
}

black_hole gen_Spin(const black_hole * const min, const black_hole * const max) {
	black_hole bh;
	bh.sp = min->sp + (max->sp - min->sp) * rand1();
	bh.th = min->th + (max->th - min->th) * rand1();
	bh.ph = min->ph + (max->ph - min->ph) * rand1();
	bh.sx = bh.sp * sqrt(1. - bh.th * bh.th) * cos(bh.ph);
	bh.sy = bh.sp * sqrt(1. - bh.th * bh.th) * sin(bh.ph);
	bh.sz = bh.sp * bh.th;
	return bh;
}

double gen_Inclination(void) {
	return min.incl + (max.incl - min.incl) * rand1();
}

double gen_Distance(void) {
	return min.dist + (max.dist - min.dist) * rand1();
}

double gen_Declination(void) {
	return min.dec + (max.dec - min.dec) * rand1();
}

double gen_Polarization(void) {
	return min.pol + (max.pol - min.pol) * rand1();
}

double gen_Phi(void) {
	return min.phi + (max.phi - min.phi) * rand1();
}

binary_system gen_Params() {
	binary_system now;
	now.bh1 = gen_Spin(&min.bh1, &max.bh1);
	now.bh2 = gen_Spin(&min.bh2, &max.bh2);
	gen_m1_m2(&now.bh1.m, &now.bh2.m, &now.M, &now.eta);
	now.incl = gen_Inclination();
	now.dist = gen_Distance();
	now.dec = gen_Declination();
	now.pol = gen_Polarization();
	now.phi = gen_Phi();
	return now;
}

void set(binary_system * const dest, const binary_system * const source) {
	dest->bh1.sp = source->bh1.sp;
	dest->bh1.sx = source->bh1.sx;
	dest->bh1.sy = source->bh1.sy;
	dest->bh1.sz = source->bh1.sz;
	dest->bh1.th = source->bh1.th;
	dest->bh1.ph = source->bh1.ph;
	dest->bh1.m = source->bh1.m;
	dest->bh2.sp = source->bh2.sp;
	dest->bh2.sx = source->bh2.sx;
	dest->bh2.sy = source->bh2.sy;
	dest->bh2.sz = source->bh2.sz;
	dest->bh2.th = source->bh2.th;
	dest->bh2.ph = source->bh2.ph;
	dest->bh2.m = source->bh2.m;
	dest->incl = source->incl;
	dest->dist = source->dist;
	dest->eta = source->eta;
	dest->M = source->M;
	dest->dec = source->dec;
	dest->phi = source->phi;
	dest->pol = source->pol;
}
