/**
 * @file util_math.h
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <time.h>
#include <math.h>
#include "util.h"

#define SQR(A) ((A)*(A))///<a
typedef struct {
	const double MINUTE_TO_SECOND;
	const double SECOND_TO_MINUTE;
	const double HOUR_TO_MINUTE;
	const double MINUTE_TO_HOUR;
	const double HOUR_TO_SECOND;
	const double SECOND_TO_HOUR;
	const double DAY_TO_HOUR;
	const double HOUR_TO_DAY;
	const double DAY_TO_MINUTE;
	const double MINUTE_TO_DAY;
	const double DAY_TO_SECOND;
	const double SECOND_TO_DAY;
} TIME_CONVERSION_CONSTANTS;

typedef struct {
	const double RADIAN_TO_DEGREE;
	const double DEGREE_TO_RADIAN;
	const double HOUR_TO_RADIAN;
	const double RADIAN_TO_HOUR;
	const double HOUR_TO_DEGREE;
	const double DEGREE_TO_HOUR;
	const double MINUTE_TO_RADIAN;
	const double RADIAN_TO_MINUTE;
	const double SECOND_TO_RADIAN;
	const double RADIAN_TO_SECOND;
} CONVERSION_CONSTANTS;

extern const TIME_CONVERSION_CONSTANTS TIME_CONVERSION_CONSTANT;

extern const CONVERSION_CONSTANTS CONVERSION_CONSTANT;

/**  	Returns a random number between [0,1).
 * Use srand() beforhand.
 * @return the random number
 */
double rand1(void);

/**		Returns a random number between [0, n).
 * Use srand() beforhand.
 * @param[in]	n	: the upper limit of the generated numbers.
 * @return the random number.
 */
double randn(double n);

/**		Returns a random number between [lower, upper).
 * Use srand() beforhand.
 * @param[in]	lower	: the lowest possible number
 * @param[in]	upper	: the upper limit of the generated numbers.
 * @return the random number.
 */
double randnk(double lower, double upper);

/**		Returns the smallest power of two no less than num.
 * @param[in]	num	:
 * @return
 */
long ceil_po2(double num);

/**		Returns the largest power of two not greater than num.
 * @param[in]	num	:
 * @return
 */
long floor_po2(double num);

/**		Rounds to power of two.
 * @param[in]	num	:
 * @return
 */
long round_po2(double num);

double time_To_Radian(double hour, double minute, double second);

double time_To_Degree(double hour, double minute, double second);

//	vector functions

/**	Returns the dot product.
 * @param[in]	vec1 : first vector
 * @param[in]	vec2 : second vector
 * @return	the dot product
 */
double scalar_Product(double vec1[], double vec2[]);
#define SCALAR_PRODUCT(vec1,vec2) vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]///<a

/** Calculates the cross product.
 * @param[out]	vec	: the cross product
 * @param[in]	vec1: the first vector
 * @param[in]	vec2: the second vector
 */
void vector_Product(double vec[], double vec1[], double vec2[]);
#define VECTOR_PRODUCT(vec,vec1,vec2)			\
	vec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];	\
	vec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];	\
	vec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];///<a

/**	Calculates the unity vector and returns the length of the vector.
 * @param[out]	uni : the unity vector
 * @param[in]	vec : the vector
 * @return	: the length of the vector
 */
double unity_Vector(double uni[], double vec[]);

/**	Calculates the sum of two vectors
 * @param[out]	sum	: the sum of the vectors
 * @param[in]	vec1: the first vector
 * @param[in]	vec2: the second vector
 */
void add_Vector(double sum[], double vec1[], double vec2[]);
#define ADD_VECTOR(sum,vec1,vec2)	\
	sum[0] = vec1[0] + vec2[0];		\
	sum[1] = vec1[1] + vec2[1];		\
	sum[2] = vec1[2] + vec2[2];///<a

/** Returns the length of the vector.
 * @param[in] vec : the vector
 * @return	: the length of the vector
 */
double length_Vector(double vec[]);
#define LENGTH_VECTOR(vec)	sqrt(SCALAR_PRODUCT(vec, vec));///<a

/** Returns the angle between the two vectors.
 * @param[in] vec1 : the first vector
 * @param[in] vec2 : the second vector
 * @return	: the angle
 */
double angle_Vector(double vec1[], double vec2[]);
#define	ANGLE_VECTOR(vec1,vec2)	\
	acos(SCALAR_PRODUCT(vec1, vec2) / (LENGTH_VECTOR(vec1), LENGTH_VECTOR(vec2)));///<a

/** Returns the triple product of the vectors
 * @param[in] vec1 : first vector
 * @param[in] vec2 : second vector
 * @param[in] vec3 : third vector
 * @return	: the triple product
 */
double triple_Product(double vec1[], double vec2[], double vec3[]);
#define TRIPLE_PRODUCT(vec1,vec2,vec3)					\
	vec1[0]*vec2[1]*vec3[2] + vec1[1]*vec2[2]*vec3[0] + \
	vec1[2]*vec2[0]*vec3[1] - vec1[0]*vec2[2]*vec3[1] - \
	vec1[1]*vec2[0]*vec3[2] - vec1[2]*vec2[1]*vec3[0];///<a

#endif	// UTIL_MATH_H
