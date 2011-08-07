/**
 * @file detector.h
 * @author László Veréb
 * @date 2010.03.26.
 * @brief Variables and functions for detectors.
 * @see gr-qc/0008066 [1]
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include "util.h"

/** Egyéb
 */
typedef struct DETECTOR_CONSTANTS {
	const double GPS_JAN1ST2000_MIDNIGHT; ///<x
	const double LEAP_SECONDS_AT_JAN1ST2000; ///<x
	const double SECONDS_FROM_GPS_JAN1ST2000_MIDNIGHT_TO_GPS_MAY1ST2005; ///<x
	const double SECONDS_FROM_GPS_JAN1ST2000_MIDNIGHT_TO_GPS_SEP1ST2008; ///<x
	const double CENTURIES_TO_DAYS; ///<x
	const double DAYS_TO_CENTURIES; ///<x
	const double COEFFICIENT_FOR_DEGREE0; ///<x
	const double COEFFICIENT_FOR_DEGREE1; ///<x
	const double COEFFICIENT_FOR_DEGREE2; ///<x
	const double COEFFICIENT_FOR_DEGREE3; ///<x
} DETECTOR_CONSTANTS;

/**	Contains various constants used for detectors.
 */
extern const DETECTOR_CONSTANTS DETECTOR_CONSTANT;

/** The detectors
 */
typedef enum detector_Enum {
	LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT, NUMBER_OF_DETECTORS,
} detector_Enum;

/** Defines a detector
 */
typedef struct detector_table {
	detector_Enum id; ///< id of the detector
	char* name; ///< name of the detector
	double direction_Vector_X[3]; ///< x arm of the detector
	double direction_Vector_Y[3]; ///< y arm of the detector
	double location[3]; ///< location of the detector
} Detector_Table;

/**	Contains variables for the antenna beam pattern
 */
typedef struct {
	double declination; ///< in radians
	double polarization; ///< in radians
	double right_Ascention; ///< in radians
	double gmst; ///< in radians
	double greenwich_Hour_Angle; ///< in radians
	double antenna_Beam_Pattern[2]; ///< \f$F_+, F_\times\f$
} antenna_Func;

/**	Returns the id of the named detector.
 * @param[in] name	: the name of the detector
 * @return the id of the detector
 */
detector_Enum id_Of_Detector(char *name);

/**	Calculates the antenna pattern for the given detector.
 * @todo finish for all detectors
 * @param[in] id	: the detector
 * @param[out] F	: antenna function struct containing the antenna pattern
 */
void calc_Antenna_Pattern_For(detector_Enum id, antenna_Func *F);

/**	Gets the detector specific parameters.
 * @param[in] id	: the detector
 * @return	the parameters of the detector
 */
Detector_Table get_Detector_Table(detector_Enum id);

/**	Calculates the response matrix for a given detector parameters. Equation [1] (B6).
 * @param[in] detector	:  the detectors parameters
 * @param[out] response_Matrix
 */
void calc_Response_Matrix(Detector_Table detector, double response_Matrix[3][3]);

/**	Calculates the antenna pattern for a response matrix. Equations [1] (B1-5) (B7-8)
 * @param[in] response_Matrix	:
 * @param[out] antenna	:
 */
void calc_Antenna_Pattern_From_Response_Matrix(double response_Matrix[3][3], antenna_Func *antenna);

/** X
 * @param GPSsec
 * @return
 */
double convert_GMST_From_Seconds_To_Radians(double GPSsec);

/**	Prints the antenna parameters.
 * @param[in] file
 * @param[in] antenna
 * @param[in] format
 */
void print_Antenna_Function_Parameters(FILE*file, antenna_Func *antenna,
	OUTPUT_FORMAT_CONSTANTS *format);

/**	Prints the antenna parameters to plot.
 * @param[in] file
 * @param[in] antenna
 * @param[in] format
 */
void print_Antenna_Function_Parameters_To_Plot(FILE *file, antenna_Func *antenna,
	OUTPUT_FORMAT_CONSTANTS *format);

#endif /* DETECTOR_H_ */
