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
	fscanf(file, "%lg\n", &params->delta_Length);
	fscanf(file, "%s\n", params->approx[0]);
	fscanf(file, "%s\n", params->phase[0]);
	fscanf(file, "%hd\n", &params->amp_Code);
	fscanf(file, "%s\n", params->spin[0]);
	fscanf(file, "%s\n", params->approx[1]);
	fscanf(file, "%s\n", params->phase[1]);
	fscanf(file, "%hd\n", &params->amp_Code);
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
	fscanf(file, "%lg %lg\n", &parameters[0].dist, &parameters[1].dist);
	fscanf(file, "%lg %lg\n", &parameters[0].F.pol, &parameters[1].F.pol);
	parameters[0].F.pol *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.pol *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.alpha, &parameters[1].F.alpha);
	fscanf(file, "%lg %lg\n", &parameters[0].F.dec, &parameters[1].F.dec);
	parameters[0].F.dec *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	parameters[1].F.dec *= CONVERSION_CONSTANT.DEGREE_TO_RADIAN;
	fscanf(file, "%lg %lg\n", &parameters[0].F.gmst, &parameters[1].F.gmst);
	fclose(file);
}

void write_Waves_To_Files(Program_Parameters *prog, System_Parameters *parameters,
		signalStruct *sig) {
	assert(prog);
	assert(parameters);
	assert(sig);
	static short index[3] = { 0, 0, 0 };
	char temp[FILE_NAME_LENGTH];
	sprintf(temp, "%s", prog->folder);
	if (parameters->match_Minimax > parameters->critical_Match && (parameters->max_Length
			- parameters->min_Length) * parameters->time_Sampling < parameters->delta_Length) {
		sprintf(prog->folder, "%s/best", temp);
		write_Wave_To_File(prog, parameters, sig, index[0]);
		index[0]++;
		if (db > 148) {
			printf("% 11.6lg% 11.6lg% 11.6lg% 11.6lg\n", parameters->system[0].bh[0].chi_Amp,
					parameters->match_Typ, parameters->match_Best, parameters->match_Minimax);
		}
	} else if (parameters->match_Minimax > parameters->critical_Match) {
		sprintf(prog->folder, "%s/match", temp);
		write_Wave_To_File(prog, parameters, sig, index[1]);
		index[1]++;
		if (db > 148) {
			printf("% 11.6lg% 11.6lg% 11.6lg% 11.6lg\n", parameters->system[0].bh[0].chi_Amp,
					parameters->match_Typ, parameters->match_Best, parameters->match_Minimax);
		}
	}
}

void write_Wave_To_File(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		short index) {
	assert(prog);
	assert(parameters);
	assert(sig);
	char file_Name[FILE_NAME_LENGTH];
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number, prog->precision);
	sprintf(text, "%s%s%s", temp, temp, temp);
	sprintf(file_Name, "%s/data.txt", prog->folder);
	FILE *file = sfopen(file_Name, "a");
	if (file) {
		fprintf(file, text, parameters->system[0].M, parameters->system[0].eta,
				parameters->system[0].chirpM);
		fprintf(file, text, parameters->system[0].mu, parameters->system[0].bh[0].m,
				parameters->system[0].bh[1].m);
		fprintf(file, text, parameters->system[0].incl, parameters->system[0].dist,
				parameters->system[0].F.pol);
		fprintf(file, text, parameters->system[0].F.alpha, parameters->system[0].F.dec,
				parameters->system[0].F.gmst);
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
				parameters->system[0].bh[1].psi);
		fprintf(file, text, parameters->system[0].bh[0].chi[0], parameters->system[0].bh[0].chi[1],
				parameters->system[0].bh[0].chi[2]);
		fprintf(file, text, parameters->system[0].bh[1].chi[0], parameters->system[0].bh[1].chi[1],
				parameters->system[0].bh[1].chi[2]);
		fprintf(file, text, parameters->system[0].bh[0].theta, parameters->system[0].bh[0].varphi,
				parameters->system[0].bh[1].theta);
		fprintf(file, text, parameters->system[0].bh[1].varphi, parameters->system[0].coaPhase,
				parameters->system[0].coaTime);
		fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
		fprintf(file, text, parameters->system[1].M, parameters->system[1].eta,
				parameters->system[1].chirpM);
		fprintf(file, text, parameters->system[1].mu, parameters->system[1].bh[0].m,
				parameters->system[1].bh[1].m);
		fprintf(file, text, parameters->system[1].incl, parameters->system[1].dist,
				parameters->system[1].F.pol);
		fprintf(file, text, parameters->system[1].F.alpha, parameters->system[1].F.dec,
				parameters->system[1].F.gmst);
		fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
				parameters->system[1].bh[0].psi);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
				parameters->system[1].bh[1].psi);
		fprintf(file, text, parameters->system[1].bh[0].chi[0], parameters->system[1].bh[0].chi[1],
				parameters->system[1].bh[0].chi[2]);
		fprintf(file, text, parameters->system[1].bh[1].chi[0], parameters->system[1].bh[1].chi[1],
				parameters->system[1].bh[1].chi[2]);
		fprintf(file, text, parameters->system[1].bh[0].theta, parameters->system[1].bh[0].varphi,
				parameters->system[1].bh[1].theta);
		fprintf(file, text, parameters->system[1].bh[1].varphi, parameters->system[1].coaPhase,
				parameters->system[1].coaTime);
		fprintf(file, "%s %d %s ", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
		fprintf(file, text, parameters->match_Typ, parameters->match_Best,
				parameters->match_Minimax);
		fprintf(file, "%d\n", index);
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
	if (is_First) {
		sprintf(file_Name, "%s/wave%d.txtf", prog->folder, index);
	} else {
		sprintf(file_Name, "%s/wave%d.txt", prog->folder, index);
	}
	file = fopen(file_Name, "w");
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%s%s%s", temp, temp, temp);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[0].M, parameters->system[0].bh[0].m
			/ parameters->system[0].bh[1].m, parameters->system[0].eta);
	fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
			parameters->system[0].bh[0].psi);
	fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
			parameters->system[0].bh[1].psi);
	fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[1].M, parameters->system[1].bh[0].m
			/ parameters->system[1].bh[1].m, parameters->system[1].eta);
	fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
			parameters->system[1].bh[0].psi);
	fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
			parameters->system[1].bh[1].psi);
	fprintf(file, "%s %d %s", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
	fprintf(file, text, parameters->match_Typ, parameters->match_Best, parameters->match_Minimax);
	fprintf(file, "\n");
	long i;
	for (i = 0; i < parameters->min_Length; i++) {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
				* parameters->system[0].F.F[0] + sig->signal[H1C][i] * parameters->system[0].F.F[1]);
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	if (parameters->shorter) {
		for (; i < parameters->max_Length; i++) {
			fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
					(double) i * parameters->time_Sampling);
			fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
					* parameters->system[0].F.F[0] + sig->signal[H1C][i]
					* parameters->system[0].F.F[1]);
			fprintf(file, "\n");
		}
	} else {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, "%*s\t%*s\t%*s\t", prog->width_Of_Number_To_Plot, "",
				prog->width_Of_Number_To_Plot, "", prog->width_Of_Number_To_Plot, "");
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	fclose(file);
}

/**
 *
 * @param prog
 * @param parameters
 * @param file_Name
 */
void write_Params_To_File(Program_Parameters *prog, System_Parameters *parameters, char *file_Name) {
	assert(prog);
	assert(parameters);
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number, prog->precision);
	sprintf(text, "%s%s%s", temp, temp, temp);
	FILE *file = sfopen(file_Name, "w");
	if (file) {
		fprintf(file, text, parameters->system[0].M, parameters->system[0].eta,
				parameters->system[0].chirpM);
		fprintf(file, text, parameters->system[0].mu, parameters->system[0].bh[0].m,
				parameters->system[0].bh[1].m);
		fprintf(file, text, parameters->system[0].incl, parameters->system[0].dist,
				parameters->system[0].F.pol);
		fprintf(file, text, parameters->system[0].F.alpha, parameters->system[0].F.dec,
				parameters->system[0].F.gmst);
		fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
				parameters->system[0].bh[0].psi);
		fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
				parameters->system[0].bh[1].psi);
		fprintf(file, text, parameters->system[0].bh[0].chi[0], parameters->system[0].bh[0].chi[1],
				parameters->system[0].bh[0].chi[2]);
		fprintf(file, text, parameters->system[0].bh[1].chi[0], parameters->system[0].bh[1].chi[1],
				parameters->system[0].bh[1].chi[2]);
		fprintf(file, text, parameters->system[0].bh[0].theta, parameters->system[0].bh[0].varphi,
				parameters->system[0].bh[1].theta);
		fprintf(file, text, parameters->system[0].bh[1].varphi, parameters->system[0].coaPhase,
				parameters->system[0].coaTime);
		fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
		fprintf(file, text, parameters->system[1].M, parameters->system[1].eta,
				parameters->system[1].chirpM);
		fprintf(file, text, parameters->system[1].mu, parameters->system[1].bh[0].m,
				parameters->system[1].bh[1].m);
		fprintf(file, text, parameters->system[1].incl, parameters->system[1].dist,
				parameters->system[1].F.pol);
		fprintf(file, text, parameters->system[1].F.alpha, parameters->system[1].F.dec,
				parameters->system[1].F.gmst);
		fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
				parameters->system[1].bh[0].psi);
		fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
				parameters->system[1].bh[1].psi);
		fprintf(file, text, parameters->system[1].bh[0].chi[0], parameters->system[1].bh[0].chi[1],
				parameters->system[1].bh[0].chi[2]);
		fprintf(file, text, parameters->system[1].bh[1].chi[0], parameters->system[1].bh[1].chi[1],
				parameters->system[1].bh[1].chi[2]);
		fprintf(file, text, parameters->system[1].bh[0].theta, parameters->system[1].bh[0].varphi,
				parameters->system[1].bh[1].theta);
		fprintf(file, text, parameters->system[1].bh[1].varphi, parameters->system[1].coaPhase,
				parameters->system[1].coaTime);
		fprintf(file, "%s %hd %s ", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
		fprintf(file, text, parameters->match_Typ, parameters->match_Best,
				parameters->match_Minimax);
		fclose(file);
	} else {
		printf("Can not open file: %s, terminating!!!\n", file_Name);
		fflush(stdout);
		abort();
	}
}

void write_Waves(Program_Parameters *prog, System_Parameters *parameters, signalStruct *sig,
		char *file_Name) {
	FILE *file = fopen(file_Name, "w");
	static char temp[FILE_NAME_LENGTH];
	static char text[FILE_NAME_LENGTH];
	sprintf(temp, "%%%d.%dlg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot);
	sprintf(text, "%s%s%s", temp, temp, temp);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[0].M, parameters->system[0].bh[0].m
			/ parameters->system[0].bh[1].m, parameters->system[0].eta);
	fprintf(file, text, parameters->system[0].bh[0].chi_Amp, parameters->system[0].bh[0].kappa,
			parameters->system[0].bh[0].psi);
	fprintf(file, text, parameters->system[0].bh[1].chi_Amp, parameters->system[0].bh[1].kappa,
			parameters->system[0].bh[1].psi);
	fprintf(file, "%s %d %s\n", parameters->phase[0], parameters->amp_Code, parameters->spin[0]);
	fprintf(file, "#");
	fprintf(file, text, parameters->system[1].M, parameters->system[1].bh[0].m
			/ parameters->system[1].bh[1].m, parameters->system[1].eta);
	fprintf(file, text, parameters->system[1].bh[0].chi_Amp, parameters->system[1].bh[0].kappa,
			parameters->system[1].bh[0].psi);
	fprintf(file, text, parameters->system[1].bh[1].chi_Amp, parameters->system[1].bh[1].kappa,
			parameters->system[1].bh[1].psi);
	fprintf(file, "%s %d %s", parameters->phase[1], parameters->amp_Code, parameters->spin[1]);
	fprintf(file, text, parameters->match_Typ, parameters->match_Best, parameters->match_Minimax);
	fprintf(file, "\n");
	long i;
	for (i = 0; i < parameters->min_Length; i++) {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
				* parameters->system[0].F.F[0] + sig->signal[H1C][i] * parameters->system[0].F.F[1]);
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	if (parameters->shorter) {
		for (; i < parameters->max_Length; i++) {
			fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
					(double) i * parameters->time_Sampling);
			fprintf(file, text, sig->signal[H1P][i], sig->signal[H1C][i], sig->signal[H1P][i]
					* parameters->system[0].F.F[0] + sig->signal[H1C][i]
					* parameters->system[0].F.F[1]);
			fprintf(file, "\n");
		}
	} else {
		fprintf(file, "%*.*lg\t", prog->width_Of_Number_To_Plot, prog->precision_To_Plot,
				(double) i * parameters->time_Sampling);
		fprintf(file, "%*s\t%*s\t%*s\t", prog->width_Of_Number_To_Plot, "",
				prog->width_Of_Number_To_Plot, "", prog->width_Of_Number_To_Plot, "");
		fprintf(file, text, sig->signal[H2P][i], sig->signal[H2C][i], sig->signal[H2P][i]
				* parameters->system[0].F.F[0] + sig->signal[H2C][i] * parameters->system[0].F.F[1]);
		fprintf(file, "\n");
	}
	fclose(file);
}
