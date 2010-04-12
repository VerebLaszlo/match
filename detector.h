/**
 * @file detector.h
 * @author László Veréb
 * @date 2010.03.26.
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include "util.h"

/**
 *		Enumeration of the detectors
 */
typedef enum detector_Enum {
	LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT
} detector;

/**
 *		The function converts the detector-string to detector-enum value.
 * @param[in]	: detector string
 * @return	the detectors enum-value
 */
detector con_Det_Str_Enum(char *str);

/**
 *		The function calculates the response-matrix of the detector given by teh nx, ny vectors.
 * @param[in]	nx : unity-vector of the detectors x-arm
 * @param[in]	ny : unity-vector of the detectors y-arm
 * @param[out]	rm : the detectors response-matri
 */
void calc_Response_Matrix(const double nx[3], const double ny[3], double rm[3][3]);

/**
 *		The function calulates the respons-function of the detector given by the response-matrix.
 * @param[in]	D	: the detectors response-matrix
 * @param[in]	dec	: the declination of the source
 * @param[in]	phi	: \todo ez mi
 * @param[in]	pol	: the sources polarization angle
 * @param[out]	fp	: F+
 * @param[out]	fc	: Fx
 */
void calc_Response(double D[3][3], double dec, double phi, double pol, double *fp, double *fc);

/**
 *		The function calulates the given detectors response-function. Currentli supports only three.
 *	\bug	Még be kell fejezni a többi detektorra
 * @param[in]	det	: the detectors ID from the enumeration, (LL : Livingston, LH : Hanford, VIRGO: Virgo)
 * @param[in]	dec	: the declination of the source
 * @param[in]	phi	: \todo ez mi
 * @param[in]	pol	: the sources polarization angle
 * @param[out] fp	: F+
 * @param[out] fc	: Fx
 */
void calc_Response_For_Detector(detector det, double dec, double phi, double pol, double *fp, double *fc);

#endif /* DETECTOR_H_ */
