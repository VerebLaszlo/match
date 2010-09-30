/**
 * @file util_math.c
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#include "util_math.h"

double rand1(void) {
	return (double)rand() / ((double)RAND_MAX + 1.);
}

double randn(double n) {
	return n * (double)rand() / ((double)RAND_MAX + 1.);
}

double randnk(double lower, double upper) {
	return (lower + upper) * (double)rand() / ((double)RAND_MAX + 1.) - lower;
}

long ceil_po2(double num) {
	register double temp = log(num) / M_LN2;
#ifdef __USE_ISOC99
	return (long) exp2(ceil(temp));
#endif
	//return (long) pow(2., ceil(temp));
	return (long) exp(ceil(temp) * M_LN2);
}

//	vector functions

double scalar_Product(double vec1[], double vec2[]) {
	return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

void vector_Product(double vec[], double vec1[], double vec2[]) {
	vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
	vec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
	vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

double unity_Vector(double uni[], double vec[]) {
	double length = length_Vector(vec);
	uni[0] = vec[0] / length;
	uni[1] = vec[1] / length;
	uni[2] = vec[2] / length;
	return length;
}

void add_Vector(double sum[], double vec1[], double vec2[]) {
	sum[0] = vec1[0] + vec2[0];
	sum[1] = vec1[1] + vec2[1];
	sum[2] = vec1[2] + vec2[2];
}

double length_Vector(double vec[]) {
	return sqrt(scalar_Product(vec, vec));
}

double angle_Vector(double vec1[], double vec2[]) {
	return acos(scalar_Product(vec1, vec2) /
			(length_Vector(vec1), length_Vector(vec2)));
}

double triple_Product(double vec1[], double vec2[], double vec3[]) {
	return	vec1[0]*vec2[1]*vec3[2] + vec1[1]*vec2[2]*vec3[0] +
			vec1[2]*vec2[0]*vec3[1] - vec1[0]*vec2[2]*vec3[1] -
			vec1[1]*vec2[0]*vec3[2] - vec1[2]*vec2[1]*vec3[0];
}

