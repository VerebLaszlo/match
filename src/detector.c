/**
 * @file detector.c
 * @author László Veréb
 * @date 2011.08.16.
 */

#include <assert.h>
#include <stdlib.h>
#include "detector.h"
#include "test.h"
#include "util_math.h"

const DetectorTable detectors[] = { //
	{ LL, "LL", //
		{ //
		{ -0.9546, -0.1416, -0.2622 }, //
		{ +0.2977, -0.4879, -0.8205 } //
		}, //
		{ -0.074276, -5.496284, 3.224257 }, //
	},//
	{ LH, "LH", //
		{ //
		{ -0.2239, 0.7998, 0.5569 }, //
		{ -0.9140, 0.0261, -0.4049 }, //
		},//
		{ -2.161415, -3.834695, 4.600350 }, //
	},//
	{ VIRGO, "VIRGO", //
		{ //
		{ -0.7005, 0.2085, 0.6826 }, //
		{ -0.0538, -0.9691, 0.2408 }, //
		},//
		{ 4.546374, 0.842990, 4.378577 }, //
	},//
	{ GEO600, "GEO600", //
		{ //
		{ -0.6261, -0.5522, 0.5506 }, //
		{ -0.4453, 0.8665, 0.2255 }, //
		},//
		{ 3.856310, 0.666599, 5.019641 }, //
	},//
	{ TAMA20, "TAMA20", //
		{ //
		{ 0.7727, 0.2704, 0.5744 }, //
		{ -0.1451, -0.8056, 0.5744 }, //
		},//
		{ -3.946416, 3.365795, 3.699409 }, //
	},//
	{ TAMA300, "TAMA300", //
		{ //
		{ 0.6490, 0.7608, 0.0000 }, //
		{ -0.4437, 0.3785, -0.8123 }, //
		},//
		{ -3.946409, 3.366259, 3.699151 }, //
	},//
	{ GLASGOW, "GLASGOW", //
		{ //
		{ -0.4534, -0.8515, 0.2634 }, //
		{ 0.6938, -0.5227, -0.4954 }, //
		},//
		{ 3.576830, -0.267688, 5.256335 }, //
	},//
	{ ISAS100, "ISAS100", //
		{ //
		{ 0.7634, 0.2277, 0.6045 }, //
		{ 0.1469, 0.8047, -0.5752 }, //
		},//
		{ -3.947704, 3.375234, 3.689488 }, //
	},//
	{ MPQ, "MPQ", //
		{ //
		{ -0.7304, 0.3749, 0.5709 }, //
		{ 0.2027, 0.9172, -0.3430 }, //
		},//
		{ 4.167725, 0.861577, 4.734691 }, //
	},//
	{ CIT, "CIT", //
		{ //
		{ -0.2648, -0.4953, -0.8274 }, //
		{ 0.8819, -0.4715, 0.0000 }, //
		},//
		{ -2.490650, -4.658700, 3.562064 }, //
	}//
	};

/**	Returns the id of the named detector.
 * @param[in] name : the name of the detector
 * @return the id of the detector
 */
static DetectorID getDetectorIDOf(const char *name) {
	BACKUP_DEFINITION_LINE();
	for (short i = 0; i < NUMBER_OF_DETECTORS; i++) {
		if (!strcmp(name, detectors[i].name)) {
			SAVE_FUNCTION_FOR_TESTING();
			return detectors[i].id;
		}
	}
	fprintf(stderr, "There is not such detector as: %s!", name);
	exit(EXIT_FAILURE);
}

/**	Gets the detector specific parameters.
 * @param[in] id : the detector
 * @return	the parameters of the detector
 */
static DetectorTable getDetectorTable(DetectorID id) {
	BACKUP_DEFINITION_LINE();
	for (short i = 0; i < NUMBER_OF_DETECTORS; i++) {
		if (detectors[i].id == id) {
			SAVE_FUNCTION_FOR_TESTING();
			return detectors[i];
		}
	} //
	SAVE_FUNCTION_FOR_TESTING();
	return detectors[0];
}

/**	Calculates the response matrix for a given detector parameters. Equation [1] (B6).
 * @param[out] responseMatrix : the response matrix of the detector
 * @param[in]  detector		  : the detectors parameters
 */
static void calcResponseMatrix(double responseMatrix[DIMENSION][DIMENSION], DetectorTable detector) {
	BACKUP_DEFINITION_LINE();
	for (short i = 0; i < DIMENSION; i++) {
		responseMatrix[i][i] = (detector.arm[X][i] * detector.arm[X][i]
			- detector.arm[Y][i] * detector.arm[Y][i]) / 2.0;
		for (short j = i + 1; j < DIMENSION; j++) {
			responseMatrix[i][j] = responseMatrix[j][i] = (detector.arm[X][i] * detector.arm[X][j]
				- detector.arm[Y][i] * detector.arm[Y][j]) / 2.0;
		}
	} //
	SAVE_FUNCTION_FOR_TESTING();
}

/**	Calculates the antenna pattern for a response matrix. Equations [1] (B1-5) (B7-8)
 * @param[out] antenna		   :
 * @param[in]  responseMatrix :
 */
static void calcAntennaPatternFromResponseMatrix(DetectorParameters *antenna,
	double responseMatrix[DIMENSION][DIMENSION]) {
	BACKUP_DEFINITION_LINE();
	double VEC1[DIMENSION];
	double VEC2[DIMENSION];
	const double cosGHA = cosGood(antenna->greenwichHourAngle);
	const double sinGHA = sinGood(antenna->greenwichHourAngle);
	const double cosDeclination = cosGood(antenna->declination);
	const double sinDeclination = sinGood(antenna->declination);
	const double cosPolarization = cosGood(antenna->polarization);
	const double sinPolarization = sinGood(antenna->polarization);
	VEC1[X] = -cosPolarization * sinGHA - sinPolarization * cosGHA * sinDeclination;
	VEC1[Y] = -cosPolarization * cosGHA + sinPolarization * sinGHA * sinDeclination;
	VEC1[Z] = sinPolarization * cosDeclination;
	VEC2[X] = sinPolarization * sinGHA - cosPolarization * cosGHA * sinDeclination;
	VEC2[Y] = sinPolarization * cosGHA + cosPolarization * sinGHA * sinDeclination;
	VEC2[Z] = cosPolarization * cosDeclination;
	antenna->antennaBeamPattern[0] = antenna->antennaBeamPattern[1] = 0.0;
	for (short dim = 0; dim < DIMENSION; dim++) {
		const double DX = responseMatrix[dim][X] * VEC1[X] + responseMatrix[dim][Y] * VEC1[Y]
			+ responseMatrix[dim][Z] * VEC1[Z];
		const double DY = responseMatrix[dim][X] * VEC2[X] + responseMatrix[dim][Y] * VEC2[Y]
			+ responseMatrix[dim][Z] * VEC2[Z];
		antenna->antennaBeamPattern[0] += VEC1[dim] * DX - VEC2[dim] * DY;
		antenna->antennaBeamPattern[1] += VEC1[dim] * DY + VEC2[dim] * DX;
	} //
	SAVE_FUNCTION_FOR_TESTING();
}

void calcAntennaPatternFor(DetectorID id, DetectorParameters *parameter) {
	BACKUP_DEFINITION_LINE();
	double response_Matrix[3][3];
	parameter->greenwichHourAngle = parameter->greenwichMeanSiderealTime
		- parameter->rightAscention;
	DetectorTable detector = getDetectorTable(id);
	calcResponseMatrix(response_Matrix, detector);
	calcAntennaPatternFromResponseMatrix(parameter, response_Matrix);
	SAVE_FUNCTION_FOR_TESTING();
}

void generateDetectorParameters(DetectorParameters *detector, DetectorParameters limits[]) {
	BACKUP_DEFINITION_LINE(); //
	assert(detector);
	assert(limits);
	detector->declination = randomBetween(limits[MIN].declination, limits[MAX].declination);
	detector->polarization = randomBetween(limits[MIN].polarization, limits[MAX].polarization);
	detector->rightAscention = randomBetween(limits[MIN].rightAscention,
		limits[MAX].rightAscention);
	detector->greenwichMeanSiderealTime = randomBetween(limits[MIN].greenwichMeanSiderealTime,
		limits[MAX].greenwichMeanSiderealTime);
	SAVE_FUNCTION_FOR_TESTING();
}

void printDetectorParameters(FILE *file, DetectorParameters *detector, OutputFormat *format) {
	BACKUP_DEFINITION_LINE();
	ushort number = 4;
	ushort length = number * format->widthWithSeparator;
	char formatString[length];
	setFormat(formatString, number, format);
	fprintf(file, formatString, detector->declination, detector->polarization,
		detector->rightAscention, detector->greenwichMeanSiderealTime);
	number = 3;
	length = number * format->widthWithSeparator;
	setFormatEnd(formatString, number, format);
	fprintf(file, formatString, detector->greenwichHourAngle, detector->antennaBeamPattern[0],
		detector->antennaBeamPattern[1]);
	SAVE_FUNCTION_FOR_TESTING();
}

#ifdef TEST

static bool isOK_getDetectorIDOf(void) {
	DetectorID id[] = {LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT,
		NUMBER_OF_DETECTORS,};
	const char *name[] = {"LL", "LH", "VIRGO", "GEO600", "TAMA20", "TAMA300", "GLASGOW", "ISAS100",
		"MPQ", "CIT",};
	ushort detector = 0;
	while (id[detector] != NUMBER_OF_DETECTORS) {
		SAVE_FUNCTION_CALLER();
		if (id[detector] != getDetectorIDOf(name[detector])) {
			PRINT_ERROR();
			return false;
		}
		detector++;
	}
	detector = 0;
	while (id[detector] != NUMBER_OF_DETECTORS) {
		SAVE_FUNCTION_CALLER();
		if (id[detector + 1] == getDetectorIDOf(name[detector])) {
			PRINT_ERROR();
			return false;
		}
		detector++;
	} //
	PRINT_OK();
	return true;
}

static bool isOK_getDetectorTable(void) {
	DetectorID id[] = {LL, LH, VIRGO, GEO600, TAMA20, TAMA300, GLASGOW, ISAS100, MPQ, CIT,
		NUMBER_OF_DETECTORS,};
	const char *name[] = {"LL", "LH", "VIRGO", "GEO600", "TAMA20", "TAMA300", "GLASGOW", "ISAS100",
		"MPQ", "CIT",};
	ushort detector = 0;
	DetectorTable table;
	while (id[detector] != NUMBER_OF_DETECTORS) {
		SAVE_FUNCTION_CALLER();
		table = getDetectorTable(id[detector]);
		if (strcmp(table.name, name[detector])) {
			PRINT_ERROR();
			return false;
		}
		detector++;
	}
	detector = 0;
	while (id[detector] != NUMBER_OF_DETECTORS) {
		SAVE_FUNCTION_CALLER();
		table = getDetectorTable(id[detector + 1]);
		if (!strcmp(table.name, name[detector])) {
			PRINT_ERROR();
			return false;
		}
		detector++;
	} //
	PRINT_OK();
	return true;
}

static bool isOK_calcResponseMatrix(void) {
	ushort detector = LL;
	double matrix[DIMENSION][DIMENSION];
	double result[GEO600][DIMENSION][DIMENSION] = { //
		{ //
			{	+0.4113180, +0.1402100, +0.2472790}, //
			{	+0.1402100, -0.1089980, -0.1815970}, //
			{	+0.2472790, -0.1815970, -0.3022360}, //
		},
		{ //
			{	-0.3926320, -0.0776099, -0.2473840}, //
			{	-0.0776099, +0.3194990, +0.2279880}, //
			{	-0.2473840, +0.2279880, +0.0730968}, //
		},
		{ //
			{	+0.2439030, -0.0990959, -0.2326030}, //
			{	-0.0990959, -0.4478410, +0.1878410}, //
			{	-0.2326030, +0.1878410, +0.2039790}, //
		}, //
	};
	while (detector < GEO600) {
		SAVE_FUNCTION_CALLER();
		calcResponseMatrix(matrix, detectors[detector]);
		for (ushort dim1 = X; dim1 < DIMENSION; dim1++) {
			for (ushort dim2 = X; dim2 < DIMENSION; dim2++) {
				if (!isNear(matrix[dim1][dim2], result[detector][dim1][dim2], 4.849999999946e-07)) {
					PRINT_ERROR();
					return false;
				}
			}
		}
		detector++;
	} //
	PRINT_OK();
	return true;
}

static bool isOK_calcAntennaPatternFromResponseMatrix(void) {
	if (!isOK_calcResponseMatrix()) {
		return false;
	}
	double matrix[DIMENSION][DIMENSION];
	DetectorParameters antenna;
	antenna.declination = antenna.polarization = antenna.greenwichHourAngle = 1.0;
	calcResponseMatrix(matrix, detectors[LL]);
	double epsilon = 3.406948470563443e-07;
	SAVE_FUNCTION_CALLER();
	calcAntennaPatternFromResponseMatrix(&antenna, matrix);
	if (!isNear(antenna.antennaBeamPattern[0], -0.1663658543172435722, epsilon)
		|| !isNear(antenna.antennaBeamPattern[1], -0.79869393295667401311, epsilon)) {
		PRINT_ERROR();
		return false;
	} //
	PRINT_OK();
	return true;
}

static bool isOK_calcAntennaPatternFor(void) {
	DetectorParameters params;
	params.greenwichMeanSiderealTime = 2.0;
	params.rightAscention = 1.0;
	params.declination = 1.0;
	params.polarization = 1.0;
	double epsilon = 3.406948470563443e-07;
	SAVE_FUNCTION_CALLER();
	calcAntennaPatternFor(LL, &params);
	if (!isNear(params.antennaBeamPattern[0], -0.1663658543172435722, epsilon)
		|| !isNear(params.antennaBeamPattern[1], -0.79869393295667401311, epsilon)) {
		PRINT_ERROR();
		return false;
	} //
	PRINT_OK();
	return true;
}

static bool isOK_generateDetectorParameters(void) {
	if (!areUtilMathFunctionsOK()) {
		return false;
	}
	DetectorParameters params, limits[2];
	limits[MAX].polarization = 10.0 + (limits[MIN].polarization = 0.0);
	limits[MAX].declination = 10.0 + (limits[MIN].declination = 0.0);
	limits[MAX].rightAscention = 10.0 + (limits[MIN].rightAscention = 0.0);
	limits[MAX].greenwichMeanSiderealTime = 10.0 + (limits[MIN].greenwichMeanSiderealTime = 0.0);
	SAVE_FUNCTION_CALLER();
	generateDetectorParameters(&params, limits);
	if (limits[MIN].polarization > params.polarization
		|| params.polarization > limits[MAX].polarization
		|| limits[MIN].declination > params.declination
		|| params.declination > limits[MAX].declination
		|| limits[MIN].rightAscention > params.rightAscention
		|| params.rightAscention > limits[MAX].rightAscention
		|| limits[MIN].greenwichMeanSiderealTime > params.greenwichMeanSiderealTime
		|| params.greenwichMeanSiderealTime > limits[MAX].greenwichMeanSiderealTime) {
		PRINT_ERROR();
		return false;
	} //
	PRINT_OK();
	return true;
}

bool areDetectorFunctionsGood(void) {
	bool isOK = true;
	if (!isOK_getDetectorIDOf()) {
		isOK = false;
	} else if (!isOK_getDetectorTable()) {
		isOK = false;
	} else if (!isOK_calcAntennaPatternFromResponseMatrix()) {
		isOK = false;
	} else if (!isOK_calcAntennaPatternFor()) {
		isOK = false;
	} else if (!isOK_generateDetectorParameters()) {
		isOK = false;
	}
	if (isOK) {
		PRINT_OK_FILE();
	} else {
		PRINT_ERROR_FILE();
	}
	return isOK;
}

#endif	// TEST
/*
 const DETECTOR_CONSTANTS DETECTOR_CONSTANT = { 630720013.0, 32.0, 189388800.0, 284083201.0, 36525.0,
 1.0 / 36525.0, 24110.54841, 8640184.812866,
 0.093104, 6.2e-6, };
 */

/*
 double convert_GMST_From_Seconds_To_Radians(double GPSsec) {
 //double GPS_Jan1st2000midnight = 630720013.0;
 double leapseconds = DETECTOR_CONSTANT.LEAP_SECONDS_AT_JAN1ST2000; // at Jan 1st 2000
 double days_Since_Jan1st2000, centuries, seconds_In_Current_Day, mean_Sidereal_Time;
 if (GPSsec
 > (DETECTOR_CONSTANT.GPS_JAN1ST2000_MIDNIGHT
 + DETECTOR_CONSTANT.SECONDS_FROM_GPS_JAN1ST2000_MIDNIGHT_TO_GPS_MAY1ST2005)) {
 leapseconds += 1.0;
 }
 if (GPSsec
 > (DETECTOR_CONSTANT.GPS_JAN1ST2000_MIDNIGHT
 + DETECTOR_CONSTANT.SECONDS_FROM_GPS_JAN1ST2000_MIDNIGHT_TO_GPS_SEP1ST2008)) {
 leapseconds += 1.0;
 }
 if (GPSsec < DETECTOR_CONSTANT.GPS_JAN1ST2000_MIDNIGHT) {
 printf(" : WARNING: GMST's before 1.1.2000 may be inaccurate! \n");
 printf(" :          (requested: GMST(GPS=%.3fs))\n", GPSsec);
 }
 // time since Jan 1st 2000 (0:00h)
 double seconds_Since_Jan1st2000 = (GPSsec - DETECTOR_CONSTANT.GPS_JAN1ST2000_MIDNIGHT)
 + (leapseconds - DETECTOR_CONSTANT.LEAP_SECONDS_AT_JAN1ST2000);
 days_Since_Jan1st2000 = floor(seconds_Since_Jan1st2000 * TIME_CONVERSION_CONSTANT.SECOND_TO_DAY)
 - 0.5;
 seconds_In_Current_Day = fmod(seconds_Since_Jan1st2000, TIME_CONVERSION_CONSTANT.DAY_TO_SECOND);
 centuries = days_Since_Jan1st2000 * DETECTOR_CONSTANT.DAYS_TO_CENTURIES;
 mean_Sidereal_Time = DETECTOR_CONSTANT.COEFFICIENT_FOR_DEGREE0
 + (centuries
 * (DETECTOR_CONSTANT.COEFFICIENT_FOR_DEGREE1
 + centuries
 * (DETECTOR_CONSTANT.COEFFICIENT_FOR_DEGREE2
 + centuries * DETECTOR_CONSTANT.COEFFICIENT_FOR_DEGREE3)));
 double coordinated_Universal_Time = mean_Sidereal_Time
 + seconds_In_Current_Day * 1.002737909350795; // (UTC day is 1.002 * MST day)
 coordinated_Universal_Time = fmod(mean_Sidereal_Time * TIME_CONVERSION_CONSTANT.SECOND_TO_DAY, 1.0);
 coordinated_Universal_Time *= 2.0 * M_PI;
 return coordinated_Universal_Time;
 }
 */

/*
 void print_Antenna_Function_Parameters(FILE *file, Detector *antenna,
 OUTPUT_FORMAT_CONSTANTS *format) {
 char format_String[4 * format->width_Of_Number_Width_Separator];set_Format_For(format_String, 4, format);
 fprintf(file, format_String, antenna->polarization, antenna->declination,
 antenna->right_Ascention, antenna->gmst);
 }

 void print_Antenna_Function_Parameters_To_Plot(FILE*file, Detector *antenna,
 OUTPUT_FORMAT_CONSTANTS *format) {
 char format_String[4 * format->width_Of_Number_To_Plot_Width_Separator];set_Plot_Format_For(format_String, 4, format);
 fprintf(file, format_String, antenna->polarization, antenna->declination,
 antenna->right_Ascention, antenna->gmst);
 }
 */
