/**
 * @file io_handler.c
 *
 * @date Mar 2, 2011
 * @author vereb
 */

#include "io_handler.h"

void read_Program_Parameters(Program_Parameters *parameters, System_Parameters *params,
		char *file_Name) {
	assert(parameters);
	assert(params);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%ld\n", &parameters->number_Of_Runs);
	fscanf(file, "%hd\n", &parameters->precision);
	parameters->width_Of_Number = parameters->precision + EXTRA_CHARACTERS;
	fscanf(file, "%hd\n", &parameters->precision_To_Plot);
	parameters->width_Of_Number_To_Plot = parameters->precision_To_Plot + EXTRA_CHARACTERS;
	fscanf(file, "%s\n", parameters->folder);
	fscanf(file, "%lg\n", &params->min_Match);
	fscanf(file, "%lg\n", &params->max_Spin);
	fscanf(file, "%lg\n", &params->spin_Step);
	fscanf(file, "%lg\n", &params->freq_Sampling);
	params->time_Sampling = 1. / params->freq_Sampling;
	fscanf(file, "%lg\n", &params->freq_Initial);
	fscanf(file, "%lg\n", &params->freq_Max);
	fscanf(file, "%lg\n", &params->delta_Length);
	fscanf(file, "%s\n", params->approx[0]);
	fscanf(file, "%s\n", params->phase[0]);
	fscanf(file, "%hd\n", &params->amp_Code[0]);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%hd\n", &params->amp_Code[1]);
	fscanf(file, "%s\n", params->spin[1]);
	fclose(file);
}

void read_Parameters(binary_System *parameters, char *file_Name) {
	assert(parameters);
	assert(file_Name);
	FILE *file = fopen(file_Name, "r");
	fscanf(file, "%lg %lg\n", &parameters[0].M, &parameters[1].M);
	fscanf(file, "%lg %lg\n", &parameters[0].eta, &parameters[1].eta);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].m, &parameters[1].bh[0].m);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].m, &parameters[1].bh[1].m);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].chi_Amp, &parameters[1].bh[0].chi_Amp);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].chi_Amp, &parameters[1].bh[1].chi_Amp);
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].kappa, &parameters[1].bh[0].kappa);
	parameters[0].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].kappa, &parameters[1].bh[1].kappa);
	parameters[0].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[0].psi, &parameters[1].bh[0].psi);
	parameters[0].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].bh[1].psi, &parameters[1].bh[1].psi);
	parameters[0].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].incl, &parameters[1].incl);
	parameters[0].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].dist, &parameters[1].dist);
	fscanf(file, "%lg %lg\n", &parameters[0].F.polarization, &parameters[1].F.polarization);
	parameters[0].F.polarization *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.polarization *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.right_Ascention, &parameters[1].F.right_Ascention);
	parameters[0].F.right_Ascention *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	parameters[1].F.right_Ascention *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.declination, &parameters[1].F.declination);
	parameters[0].F.declination *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.declination *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.gmst, &parameters[1].F.gmst);
	parameters[0].F.declination *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	parameters[1].F.declination *= CONVERSION_CONSTANT.SECOND_TO_RADIAN;
	fclose(file);
}

void readExactParameters(FILE *file, System_Parameters *params) {
	char name[100];
	for (short i = 0; i < 2; i++) {
		fscanf(file, "%s ", name);
		fscanf(file, "%lg ", &params->system[i].bh[0].m);
		fscanf(file, "%lg ", &params->system[i].bh[1].m);
		fscanf(file, "%lg ", &params->system[i].bh[0].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[0].kappa);
		params->system[i].bh[0].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[0].psi);
		params->system[i].bh[0].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].chi_Amp);
		fscanf(file, "%lg ", &params->system[i].bh[1].kappa);
		params->system[i].bh[1].kappa *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].bh[1].psi);
		params->system[i].bh[1].psi *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->system[i].dist);
		fscanf(file, "%lg ", &params->system[i].incl);
		params->system[i].incl *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
		fscanf(file, "%lg ", &params->freq_Sampling);
		params->time_Sampling = 1.0 / params->freq_Sampling;
		fscanf(file, "%lg ", &params->freq_Initial);
		fscanf(file, "%lg ", &params->freq_Max);
		fscanf(file, "%s ", params->phase[i]);
		fscanf(file, "%s ", params->spin[i]);
		fscanf(file, "%hd ", &params->amp_Code[i]);
		fscanf(file, "%s\n", params->approx[i]);
		params->system[i].F.declination = params->system[i].F.gmst
				= params->system[i].F.greenwich_Hour_Angle = params->system[i].F.polarization
						= params->system[i].F.right_Ascention = 0.;
		convert_Masses(&params->system[i], FROM_M1M2);
		convert_Spins(&params->system[i], FROM_KAPPA_PSI);
		calc_Antenna_Pattern_For(LH, &params->system[i].F);
	}
}
