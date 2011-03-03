/**
 * @file util_math.h
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <math.h>
#include "util.h"

/**		Calculates the square of the parameter.
 * The argument evaluates only once, so there is no undefined side effect
 * @todo examine the side effects
 * @param[in] number
 */
#define SQR(number) ({ typeof(number) _number = (number); _number*_number; })

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

/**  	Returns a pseudo-random number in the range [0,1). Use srand() beforehand to initialize the
 * random number generator.
 * @return a random number
 */
double random_From_Zero_To_One(void);

/**  	Returns a pseudo-random number in the range [0,top). Use srand() beforehand to initialize
 * the random number generator.
 * @param[in]	top	: the upper limit of the generated number.
 * @return a random number.
 */
double random_From_Zero_To(double top);

/**  	Returns a pseudo-random number in the range [bottom,top). Use srand() beforehand to
 * initialize the random number generator. It works when bottom is greater than top, too.
 * @param[in]	bottom	: the lowest possible number
 * @param[in]	top		: the upper limit of the generated numbers.
 * @return a random number.
 */
double random_Between(double bottom, double top);

/**		Returns the smallest integral value that is not less than number and is the power of two.
 * @param[in]	number	: the number in question
 * @return
 */
long greatest_Number_That_Less_Than(double number);

/**		Returns the smallest integral value that is not less than number and is the power of two.
 * @param[in]	number	: the number in question
 * @return
 */
long least_Number_That_Greater_Than(double number);

/**		Converts time into radian.
 * @param[in] hour		: hours
 * @param[in] minute	: minutes
 * @param[in] second	: seconds
 * @return the time value in radian
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
