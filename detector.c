/**
 * @file detector.c
 * @author László Veréb
 * @date 2010.03.26.
 */

#include "detector.h"

#define DETECTOR_NUMBER 10
const detector_table Detectors[] =
	{
		{ LL, "LL",
			{ -0.9546, -0.1416, -0.2622 },
			{ 0.2977, -0.4879, -0.8205 },
			{ -0.074276, -5.496284, 3.224257 } },
		{ LH, "LH",
			{ -0.2239, 0.7998, 0.5569 },
			{ -0.9140, 0.0261, -0.4049 },
			{ -2.161415, -3.834695, 4.600350 } },
		{ VIRGO, "VIRGO",
			{ -0.7005, 0.2085, 0.6826 },
			{ -0.0538, -0.9691, 0.2408 },
			{ 4.546374, 0.842990, 4.378577 } },
		{ GEO600, "GEO600",
			{ -0.6261, -0.5522, 0.5506 },
			{ -0.4453, 0.8665, 0.2255 },
			{ 3.856310, 0.666599, 5.019641 } },
		{ TAMA20, "TAMA20",
			{ 0.7727, 0.2704, 0.5744 },
			{ -0.1451, -0.8056, 0.5744 },
			{ -3.946416, 3.365795, 3.699409 } },
		{ TAMA300, "TAMA300",
			{ 0.6490, 0.7608, 0.0000 },
			{ -0.4437, 0.3785, -0.8123 },
			{ -3.946409, 3.366259, 3.699151 } },
		{ GLASGOW, "GLASGOW",
			{ -0.4534, -0.8515, 0.2634 },
			{ 0.6938, -0.5227, -0.4954 },
			{ 3.576830, -0.267688, 5.256335 } },
		{ ISAS100, "ISAS100",
			{ 0.7634, 0.2277, 0.6045 },
			{ 0.1469, 0.8047, -0.5752 },
			{ -3.947704, 3.375234, 3.689488 } },
		{ MPQ, "MPQ",
			{ -0.7304, 0.3749, 0.5709 },
			{ 0.2027, 0.9172, -0.3430 },
			{ 4.167725, 0.861577, 4.734691 } },
		{ CIT, "CIT",
			{ -0.2648, -0.4953, -0.8274 },
			{ 0.8819, -0.4715, 0.0000 },
			{ -2.490650, -4.658700, 3.562064 } } };

detector con_Det_Str_Enum(char *str) {
	int i;

	for (i = 0; i < DETECTOR_NUMBER; i++) {
		if (strcmp(str, Detectors[i].name)) {
			return i;
		}
	}
	// szerintem valami -1 -et kéne inkább visszaadni, nem abortálni a
	// programot
	fprintf(stderr, "There is not such detector as: %s!", str);
	exit(-1);

}

detector_table GetDetectorTable(enum detector_Enum id) {
	int i;

	for (i = 0; i < DETECTOR_NUMBER; i++) {
		if (Detectors[i].id == id) {
			return Detectors[i];
		}
	}
	// default
	return Detectors[DETECTOR_NUMBER - 1];
}

// megcsinálni normálisabbra!
void calc_Response_Matrix(const double nx[], const double ny[], double rm[3][3]) {
	long i, j;
	for (i = 0; i < 3; i++) {
		rm[i][i] = (nx[i] * nx[i] - ny[i] * ny[i]) / 2.;
		for (j = i + 1; j < 3; j++) {
			rm[i][j] = rm[j][i] = (nx[i] * nx[j] - ny[i] * ny[j]) / 2.;
		}
	}
}

void calc_Response(double D[3][3], double dec, double phi, double pol,
		double *fp, double *fc) {
	int i;
	double X[3];
	double Y[3];
	const double cosphi = cos(-phi);
	const double sinphi = sin(-phi);
	const double cosdec = cos(dec);
	const double sindec = sin(dec);
	const double cospsi = cos(pol);
	const double sinpsi = sin(pol);
	X[0] = -cospsi * sinphi - sinpsi * cosphi * sindec;
	X[1] = -cospsi * cosphi + sinpsi * sinphi * sindec;
	X[2] = sinpsi * cosdec;
	Y[0] = sinpsi * sinphi - cospsi * cosphi * sindec;
	Y[1] = sinpsi * cosphi + cospsi * sinphi * sindec;
	Y[2] = cospsi * cosdec;
	*fp = *fc = 0.0;
	for (i = 0; i < 3; i++) {
		const double DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
		const double DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
		*fp += X[i] * DX - Y[i] * DY;
		*fc += X[i] * DY + Y[i] * DX;
	}
}

void calc_Response_For_Detector(detector det, double dec, double phi,
		double pol, double *fp, double *fc) {
	double rm[3][3];
	detector_table dettable;

	dettable = GetDetectorTable(det);
	calc_Response_Matrix(dettable.nx, dettable.ny, rm);
	calc_Response(rm, dec, phi, pol, fp, fc);
}
