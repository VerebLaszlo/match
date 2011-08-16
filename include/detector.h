/**
 * @file detector.h
 * @author László Veréb
 * @date 2011.08.16.
 * @brief Variables and functions for detectors.
 * @see gr-qc/0008066 [1]
 */

#ifndef DETECTOR_H_
#define DETECTOR_H_

#include "util.h"

typedef enum {
	X, Y, Z, DIMENSION,
} Coordinates;

typedef enum {
	MIN, MAX,
} Codes;

/** Detector codes
 */
typedef enum DetectorID_ {
	LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT, NUMBER_OF_DETECTORS,
} DetectorID;

/** Defines a detector
 */
typedef struct DetectorTable_ {
	DetectorID id; ///< id of the detector
	const char* name; ///< name of the detector
	double arm[2][DIMENSION];	///< coordinates of the arms of the detector
	double location[DIMENSION]; ///< location of the detector
} DetectorTable;

/**	Contains variables for the antenna beam pattern
 */
typedef struct DetectorParameters_ {
	double declination; ///< in radians
	double polarization; ///< in radians
	double rightAscention; ///< in radians
	double greenwichMeanSiderealTime; ///< in radians
	double greenwichHourAngle; ///< in radians
	double antennaBeamPattern[2]; ///< \f$F_+, F_\times\f$
} DetectorParamters;

/**	Calculates the antenna pattern for the given detector.
 * @param[out] parameter : antenna function struct containing the antenna pattern
 * @param[in]  id		 : the detector
 */
void calcAntennaPatternFor(DetectorID id, DetectorParamters *parameter);

/**	Generates the detector parameters.
 * @param[out] detector : the generated parameters
 * @param[in]  limits	: the limits of the parameters
 */
void generateDetectorParameters(DetectorParamters *detector, DetectorParamters limits[]);

/** Egyéb
 */
/*typedef struct DETECTOR_CONSTANTS {
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
 */
/**	Contains various constants used for detectors.
 */
//extern const DETECTOR_CONSTANTS DETECTOR_CONSTANT;
/** X
 * @param GPSsec
 * @return
 */
//double convert_GMST_From_Seconds_To_Radians(double GPSsec);
/**	Prints the antenna parameters.
 * @param[in] file
 * @param[in] antenna
 * @param[in] format
 */
//void print_Antenna_Function_Parameters(FILE*file, Detector *antenna,
//	OUTPUT_FORMAT_CONSTANTS *format);
/**	Prints the antenna parameters to plot.
 * @param[in] file
 * @param[in] antenna
 * @param[in] format
 */
//void print_Antenna_Function_Parameters_To_Plot(FILE *file, Detector *antenna,
//	OUTPUT_FORMAT_CONSTANTS *format);
#endif /* DETECTOR_H_ */
