/*
 * @file detector.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "detector.h"
#include "util_math.h"

const DETECTOR_CONSTANTS DETECTOR_CONSTANT = { 630720013.0, 32.0, 189388800.0, 284083201.0, 36525.0,
	1.0 / 36525.0, 24110.54841, 8640184.812866, 0.093104, 6.2e-6, };

const Detector_Table detectors[] = { //
	{ LL, "LL", //
		{ -0.9546, -0.1416, -0.2622 }, //
		{ +0.2977, -0.4879, -0.8205 }, //
		{ -0.074276, -5.496284, 3.224257 } }, //
		{ LH, "LH", //
			{ -0.2239, 0.7998, 0.5569 }, //
			{ -0.9140, 0.0261, -0.4049 }, //
			{ -2.161415, -3.834695, 4.600350 } }, //
		{ VIRGO, "VIRGO", //
			{ -0.7005, 0.2085, 0.6826 }, //
			{ -0.0538, -0.9691, 0.2408 }, //
			{ 4.546374, 0.842990, 4.378577 } }, //
		{ GEO600, "GEO600", //
			{ -0.6261, -0.5522, 0.5506 }, //
			{ -0.4453, 0.8665, 0.2255 }, //
			{ 3.856310, 0.666599, 5.019641 } }, //
		{ TAMA20, "TAMA20", //
			{ 0.7727, 0.2704, 0.5744 }, //
			{ -0.1451, -0.8056, 0.5744 }, //
			{ -3.946416, 3.365795, 3.699409 } }, //
		{ TAMA300, "TAMA300", //
			{ 0.6490, 0.7608, 0.0000 }, //
			{ -0.4437, 0.3785, -0.8123 }, //
			{ -3.946409, 3.366259, 3.699151 } }, //
		{ GLASGOW, "GLASGOW", //
			{ -0.4534, -0.8515, 0.2634 }, //
			{ 0.6938, -0.5227, -0.4954 }, //
			{ 3.576830, -0.267688, 5.256335 } }, //
		{ ISAS100, "ISAS100", //
			{ 0.7634, 0.2277, 0.6045 }, //
			{ 0.1469, 0.8047, -0.5752 }, //
			{ -3.947704, 3.375234, 3.689488 } }, //
		{ MPQ, "MPQ", //
			{ -0.7304, 0.3749, 0.5709 }, //
			{ 0.2027, 0.9172, -0.3430 }, //
			{ 4.167725, 0.861577, 4.734691 } }, //
		{ CIT, "CIT", //
			{ -0.2648, -0.4953, -0.8274 }, //
			{ 0.8819, -0.4715, 0.0000 }, //
			{ -2.490650, -4.658700, 3.562064 } } //
	};

detector_Enum id_Of_Detector(char *name) {
	for (short i = 0; i < NUMBER_OF_DETECTORS; i++) {
		if (strcmp(name, detectors[i].name)) {
			return i;
		}
	}
	fprintf(stderr, "There is not such detector as: %s!", name);
	exit(EXIT_FAILURE);
}

void calc_Antenna_Pattern_For(detector_Enum id, antenna_Func *F) {
	double response_Matrix[3][3];
	F->greenwich_Hour_Angle = F->gmst - F->right_Ascention;
	Detector_Table detector = get_Detector_Table(id);
	calc_Response_Matrix(detector, response_Matrix);
	calc_Antenna_Pattern_From_Response_Matrix(response_Matrix, F);
}

Detector_Table get_Detector_Table(detector_Enum id) {
	for (short i = 0; i < NUMBER_OF_DETECTORS; i++) {
		if (detectors[i].id == id) {
			return detectors[i];
		}
	}
	return detectors[0];
}

void calc_Response_Matrix(Detector_Table detector, double response_Matrix[3][3]) {
	for (short i = 0; i < 3; i++) {
		response_Matrix[i][i] = (detector.direction_Vector_X[i] * detector.direction_Vector_X[i]
			- detector.direction_Vector_Y[i] * detector.direction_Vector_Y[i]) / 2.;
		for (short j = i + 1; j < 3; j++) {
			response_Matrix[i][j] = response_Matrix[j][i] = (detector.direction_Vector_X[i]
				* detector.direction_Vector_X[j]
				- detector.direction_Vector_Y[i] * detector.direction_Vector_Y[j]) / 2.;
		}
	}
}

void calc_Antenna_Pattern_From_Response_Matrix(double response_Matrix[3][3], antenna_Func *antenna) {
	double X[3];
	double Y[3];
	const double cosgha = cos(antenna->greenwich_Hour_Angle);
	const double singha = sin(antenna->greenwich_Hour_Angle);
	const double cosdec = cos(antenna->declination);
	const double sindec = sin(antenna->declination);
	const double cospol = cos(antenna->polarization);
	const double sinpol = sin(antenna->polarization);
	X[0] = -cospol * singha - sinpol * cosgha * sindec;
	X[1] = -cospol * cosgha + sinpol * singha * sindec;
	X[2] = sinpol * cosdec;
	Y[0] = sinpol * singha - cospol * cosgha * sindec;
	Y[1] = sinpol * cosgha + cospol * singha * sindec;
	Y[2] = cospol * cosdec;
	antenna->antenna_Beam_Pattern[0] = antenna->antenna_Beam_Pattern[1] = 0.0;
	for (short i = 0; i < 3; i++) {
		const double DX = response_Matrix[i][0] * X[0] + response_Matrix[i][1] * X[1]
			+ response_Matrix[i][2] * X[2];
		const double DY = response_Matrix[i][0] * Y[0] + response_Matrix[i][1] * Y[1]
			+ response_Matrix[i][2] * Y[2];
		antenna->antenna_Beam_Pattern[0] += X[i] * DX - Y[i] * DY;
		antenna->antenna_Beam_Pattern[1] += X[i] * DY + Y[i] * DX;
	}
}

double convert_GMST_From_Seconds_To_Radians(double GPSsec) {
	//double GPS_Jan1st2000midnight = 630720013.0;
	double leapseconds = DETECTOR_CONSTANT.LEAP_SECONDS_AT_JAN1ST2000; /* at Jan 1st 2000 */
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
	/* time since Jan 1st 2000 (0:00h) */
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
		+ seconds_In_Current_Day * 1.002737909350795; /* (UTC day is 1.002 * MST day) */
	coordinated_Universal_Time = fmod(mean_Sidereal_Time * TIME_CONVERSION_CONSTANT.SECOND_TO_DAY,
		1.0);
	coordinated_Universal_Time *= 2.0 * M_PI;
	return coordinated_Universal_Time;
}

void print_Antenna_Function_Parameters(FILE *file, antenna_Func *antenna,
	OUTPUT_FORMAT_CONSTANTS *format) {
	char format_String[4 * format->width_Of_Number_Width_Separator];set_Format_For(format_String, 4, format);
	fprintf(file, format_String, antenna->polarization, antenna->declination,
		antenna->right_Ascention, antenna->gmst);
}

void print_Antenna_Function_Parameters_To_Plot(FILE*file, antenna_Func *antenna,
	OUTPUT_FORMAT_CONSTANTS *format) {
	char format_String[4 * format->width_Of_Number_To_Plot_Width_Separator];set_Plot_Format_For(format_String, 4, format);
	fprintf(file, format_String, antenna->polarization, antenna->declination,
		antenna->right_Ascention, antenna->gmst);
}
