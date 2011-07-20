/**
 * @file util_math_tensor.h
 *
 * @date Mar 4, 2011
 * @author vereb
 */

#ifndef UTIL_MATH_TENSOR_H_
#define UTIL_MATH_TENSOR_H_

#include <assert.h>
#include <math.h>

#define NUMBER_OF_VECTOR_DIMENSION 3

typedef struct tagVector {
	double magnitude;
	double components[NUMBER_OF_VECTOR_DIMENSION];
	double unity[NUMBER_OF_VECTOR_DIMENSION];
	double azimuth;		///< \f[[0,2\pi)\f]
	double inclination;	///< \f[[0,\pi]\f]
	double elevation;	///< \f[[-\pi/2,\pi/2]\f]
} vector;

/**	Return the angle between the two vectors assuming their sizes are 3.
 * @param[in] left	: the first vector
 * @param[in] right	: the second vector
 * @return	the angle between the two vectors in radian
 */
#define	ANGLE_BETWEEN_VECTORS3(left,right) ({ \
	acos(SCALAR_PRODUCT3(left, right) / (LENGTH_OF_VECTOR3(left), LENGTH_OF_VECTOR3(right))); \
})

/**	Return the angle between the two vectors assuming their sizes are 3.
 * @param[in] left	: the first vector
 * @param[in] right	: the second vector
 * @return	the angle between the two vectors in radian
 */
double angle_Between_Vectors3(double left[], double right[]);

/**	Return the angle between the two vectors assuming their sizes are equal.
 * @param[in] left	: the first vector
 * @param[in] right	: the second vector
 * @param[in] size	: the size of the vectors
 * @return	the angle between the two vectors in radian
 */
double angle_Between_Vectors(double left[], double right[], short size);

/**	Calculates the unit vector of the given vector and returns the original vector's length
 * assuming the size of the vector is 3. @callgraph @callergraph
 * @param[in]	vector		: the vector in question
 * @param[out]	unit_Vector	: the calculated unit vector
 * @return	the length of the original vector
 */
double unit_Vector3_Of(double vector[], double unit_Vector[]);

/**	Calculates the unit vector of the given vector and returns the original vector's length.
 * @callgraph @callergraph
 * @param[in]	vector		: the vector in question
 * @param[in]	size		: the size of the vector
 * @param[out]	unit_Vector	: the calculated unit vector
 * @return	the length of the original vector
 */
double unit_Vector_Of(double vector[], short size, double unit_Vector[]);

/**	Normalizes a vector assuming that its size is 3. @callgraph @callergraph
 * @param[in]	vector					: the vector in question
 * @param[in]	normalizing_Constant	: the normalizing constant
 * @param[out]	normalized				: the normalized vector
 */
void normalize_Vector3(double vector[], double normalizing_Constant, double normalized[]);

/**	Normalizes a vector. @callgraph @callergraph
 * @param[in]	vector					: the vector in question
 * @param[in]	size					: the size of the vector
 * @param[in]	normalizing_Constant	: the normalizing constant
 * @param[out]	normalized				: the normalized vector
 */
void
normalize_Vector(double vector[], short size, double normalizing_Constant, double normalized[]);

/**	Returns the length of the vector assuming that its size is 3. @callgraph @callergraph
 * @param[in] vector	: the vector in question
 * @return	the length of the vector
 */
#define LENGTH_OF_VECTOR3(vector) \
	sqrt(DOT_PRODUCT3(vector, vector));

/**	Returns the length of the vector assuming that its size is 3. @callgraph @callergraph
 * @param[in] vector	: the vector in question
 * @return	the length of the vector
 */
double length_Of_Vector3(double vector[]);

/**	Returns the length of the vector. @callgraph @callergraph
 * @param[in] vector	: the vector in question
 * @param[in] size		: the size of the vector
 * @return	the length of the vector
 */
double length_Of(double vector[], short size);

/**	Calculates the dot product of the two vectors, assuming both their sizes are 3.
 * @param[in] left	: first vector
 * @param[in] right	: second vector
 * @return	scalar product
 */
#define DOT_PRODUCT3(left, right) ({ \
	typeof(left) _left; typeof(right)_right; \
	_left[0]*_right[0] + _left[1]*_right[1] + _left[2]*_right[2]; \
})

/**	Calculates the dot product of the two vectors, assuming both their sizes are 3. @callgraph
 * @callergraph
 * @param[in] left	: first vector
 * @param[in] right	: second vector
 * @return	scalar product
 */
double dot_Product3(double left[], double right[]);

/**	Calculates the dot product of the two vectors, assuming both have the same size (as it
 * should be). @callgraph @callergraph
 * @param[in] left	: first vector
 * @param[in] right	: second vector
 * @param[in] size	: the size of both vector
 * @return	scalar product
 */
double dot_Product(double left[], double right[], short size);

/**	Calculates the cross product of the two vectors, assuming both their sizes are 3.
 * @param[in]	left	: first vector
 * @param[in]	right	: second vector
 * @param[out]	result	: vector containing the result
 */
#define CROSS_PRODUCT3(left,right,result) ({ \
	typeof(left)_left; typeof(right)_right; typeof(result) _result; \
	_result[0] = _left[1]*_right[2] - _left[2]*_right[1]; \
	_result[1] = _left[2]*_right[0] - _left[0]*_right[2]; \
	_result[2] = _left[0]*_right[1] - _left[1]*_right[0]; \
})

/**	Calculates the cross product of the two vectors, assuming both their sizes are 3.
 * @todo write a function for general case.
 * @param[in]	left	: first vector
 * @param[in]	right	: second vector
 * @param[out]	result	: vector containing the result
 */
void cross_Product3(double left[], double right[], double result[]);

/**	Adds two vectors assuming their sizes are 3.
 * @param[in]	left	: the first vector
 * @param[in]	right	: the second vector
 * @param[out]	result	: the result
 */
#define ADD_VECTORS3(left,right, result) \
for (short i = 0; i < size; i++) { \
	result[i] = left[i] + right[i]; \
}

/**	Adds two vectors assuming their sizes are 3.
 * @param[in]	left	: the first vector
 * @param[in]	right	: the second vector
 * @param[out]	result	: the result
 */
void add_Vectors3(double left[], double right[], double result[]);

/**	Adds two vectors assuming their sizes are equal.
 * @param[in]	left	: the first vector
 * @param[in]	right	: the second vector
 * @param[in]	size	: the size of the vectors
 * @param[out]	result	: the result
 */
void add_Vectors(double left[], double right[], short size, double result[]);

/** Returns the triple product of the vectors
 * @param[in] vec1 : first vector
 * @param[in] vec2 : second vector
 * @param[in] vec3 : third vector
 * @return	the triple product
 */
#define TRIPLE_PRODUCT(vec1,vec2,vec3)					\
	vec1[0]*vec2[1]*vec3[2] + vec1[1]*vec2[2]*vec3[0] + \
	vec1[2]*vec2[0]*vec3[1] - vec1[0]*vec2[2]*vec3[1] - \
	vec1[1]*vec2[0]*vec3[2] - vec1[2]*vec2[1]*vec3[0];

/** Returns the triple product of the vectors
 * @param[in] vec1 : first vector
 * @param[in] vec2 : second vector
 * @param[in] vec3 : third vector
 * @return	the triple product
 */
double triple_Product(double vec1[], double vec2[], double vec3[]);

#endif /* UTIL_MATH_TENSOR_H_ */
