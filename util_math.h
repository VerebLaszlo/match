/**
 * @file util_math.h
 * @date 2010-08-16 10.27.00
 * @author László Veréb
 * @brief Various utilities for math.
 */

#ifndef UTIL_MATH_H
#define UTIL_MATH_H

#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#ifdef FGV_MACRO
	#define scalarProduct SCALAR_PRODUCT
	#define vectorProduct VECTOR_PRODUCT
	#define addVector ADD_VECTOR
	#define lengthVector LENGTH_VECTOR
	#define angleVector ANGLE_VECTOR
	#define tripleProduct TRIPLE_PRODUCT
#else
	#define scalarProduct scalar_Product///<a
	#define vectorProduct vector_Product///<a
	#define addVector add_Vector///<a
	#define lengthVector length_Vector///<a
	#define angleVector angle_Vector///<a
	#define tripleProduct triple_Product///<a
#endif

#define SQR(A) ((A)*(A))///<a

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

double deg_To_Rad(double deg);///<a
#define DEG_TO_RAD(deg) ((deg) * M_PI / 180.)///<a

double rad_To_Deg(double rad);///<a
#define RAD_TO_DEG(rad) ((rad) * 180. / M_PI)///<a

#define HOUR_TO_DEG(H) (15. * H)///<a
#define HOUR_TO_RAD(H) (DEG_TO_RAD(HOUR_TO_DEG(H)))///<a

#define SEC_TO_DEG(S) (15. * H / 3600.)///<a
#define SEC_TO_RAD(S) (DEG_TO_RAD(SEC_TO_DEG(S)))///<a

#define TIME_TO_DEG(H,M,S) (SEC_TO_DEG(S+60.*(M+60.*H)))///<a
#define TIME_TO_RAD(H,M,S) (DEG_TO_RAD(TIME_TO_DEG(H,M,S)))///<a

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
