/**
 * @file ioHandling.c
 * @author László Veréb
 * @date 2010.01.18.
 */

#include "ioHandling.h"

extern binary_system min, max;

program_params init_program(void) {
	FILE * file = sfopen_read("program.init");
	program_params params;
	fscanf(file, "%*s %ld\n", &params.numOfGoodMatch);
	fscanf(file, "%*s %ld\n", &params.numOfBestMatch);
	fscanf(file, "%*s %lg\n", &params.deviation);
	fscanf(file, "%*s %lg\n", &params.min_time);
	fscanf(file, "%*s %lg\n", &params.max_time);
	fscanf(file, "%*s %lg\n", &params.delta_time);
	fscanf(file, "%*s %d\n", &params.precision);
	fclose(file);
	return params;
}

void init_generator(char *mode, dpc dt, dpc freq_Min, dpc freq_Max) {
	FILE * file = sfopen_read("parameters.init");
	fscanf(file, "%*s %s\n", mode);
	int temp;
	fscanf(file, "%*s %d\n", &temp); *dt = 1. / (double)temp;
	fscanf(file, "%*s %lg\n%*s %lg\n", freq_Min, freq_Max);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.M, &max.M);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.eta, &max.eta);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh1.m, &max.bh1.m);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh2.m, &max.bh2.m);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh1.sp, &max.bh1.sp);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh2.sp, &max.bh2.sp);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh1.th, &max.bh1.th);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh2.th, &max.bh2.th);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh1.ph, &max.bh1.ph);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.bh2.ph, &max.bh2.ph);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.dist, &max.dist);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.incl, &max.incl);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.dec, &max.dec);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.pol, &max.pol);
	fscanf(file, "%*s %lg\n%*s %lg\n", &min.phi, &max.phi);
	fclose(file);
}

void write_Waveh(long index, double data[], long l, double dt) {
	char file_Name[30];
	sprintf(file_Name, "%ldh.wave", index);
	FILE *file = sfopen_write(file_Name);
	long i;
	for (i = 0; i < l; i++) {
		fprintf(file, "%lg %lg\n", i * dt, data[i]);
	}
	fclose(file);
}

void write_Wavel(long index, double data[], long l, double dt) {
	char file_Name[30];
	sprintf(file_Name, "%ldl.wave", index);
	FILE *file = sfopen_write(file_Name);
	long i;
	for (i = 0; i < l; i++) {
		fprintf(file, "%lg %lg\n", i * dt, data[i]);
	}
	fclose(file);
}

void write_Wavehx(long index, double data[], long l, double dt) {
	char file_Name[30];
	sprintf(file_Name, "x%ldh.wave", index);
	FILE *file = sfopen_write(file_Name);
	long i;
	for (i = 0; i < l; i++) {
		fprintf(file, "%lg %lg\n", i * dt, data[i]);
	}
	fclose(file);
}

void write_Wavelx(long index, double data[], long l, double dt) {
	char file_Name[30];
	sprintf(file_Name, "x%ldl.wave", index);
	FILE *file = sfopen_write(file_Name);
	long i;
	for (i = 0; i < l; i++) {
		fprintf(file, "%lg %lg\n", i * dt, data[i]);
	}
	fclose(file);
}

void write_Gen_Data(ulong i, double freq_Min, double freq_Max, double match[], double
	dev[], int good[], int length, const binary_system * const params, double dt) {
	FILE * parameters = sfopen("parameters.data", "a");
	fprintf(parameters, "%8lu ", i);
	fprintf(parameters, "%13.9g %13.9g ", params->bh1.m, params->bh2.m);
	fprintf(parameters, "%13.9g %13.9g %13.9g ", params->bh1.sx, params->bh1.sy, params->bh1.sz);
	fprintf(parameters, "%13.9g %13.9g %13.9g ", params->bh2.sx, params->bh2.sy, params->bh2.sz);
	fprintf(parameters, "%13.9g %8.9g %8.9g ", params->incl, freq_Min, freq_Max);
	fprintf(parameters, "%8.9g %12.9g %12.9g", params->dist, params->M, params->eta);
	fprintf(parameters, "%12.9g %12.9g %12.9g ", params->bh1.sp, params->bh1.ph, params->bh1.th);
	fprintf(parameters, "%12.9g %12.9g %12.9g ", params->bh2.sp, params->bh2.ph, params->bh2.th);
	fprintf(parameters, "%12.9g %12.9g %12.9g %12.9g ", params->dec, params->pol, params->phi, dt);
	int j;
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%13.9g ", match[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%3d ", good[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%13.9g ", dev[j]);
	}
	fprintf(parameters, "\n");
	fclose(parameters);
}

void write_Data_to_Plot(ulong i, double match[], double dev[], int good[], int length, const binary_system * const params) {
	FILE * to_Plot = sfopen("to_Plot.data", "a");
	fprintf(to_Plot, "%8lu ", i);
	fprintf(to_Plot, "%13.9lg %13.9lg ", params->bh1.m, params->bh2.m);
	fprintf(to_Plot, "%12.9lg %12.9lg ", params->bh1.sp, params->bh2.sp);
	fprintf(to_Plot, "%12.9lg %12.9lg ", params->bh1.ph, params->bh1.th);
	fprintf(to_Plot, "%12.9lg %12.9lg\n", params->bh2.ph, params->bh2.th);
	int j;
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%13.9g ", match[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%3d ", good[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%13.9g ", dev[j]);
	}
	fprintf(to_Plot, "\n");
	fclose(to_Plot);
}
/***/

void write_Gen_Datax(ulong i, double freq_Min, double freq_Max, double match[], double
dev[], int good[], int length, const binary_system * const params, double dt) {
	FILE * parameters;
	if (!(parameters = fopen("xparameters.data", "at"))) {
		fprintf(stderr, "Can not open the file named \"parameters.data\" for append!\n");
		exit(-1);
	}
	fprintf(parameters, "%8lu ", i);
	fprintf(parameters, "%13.9g %13.9g ", params->bh1.m, params->bh2.m);
	fprintf(parameters, "%13.9g %13.9g %13.9g ", params->bh1.sx, params->bh1.sy, params->bh1.sz);
	fprintf(parameters, "%13.9g %13.9g %13.9g ", params->bh2.sx, params->bh2.sy, params->bh2.sz);
	fprintf(parameters, "%13.9g %8.9g %8.9g ", params->incl, freq_Min, freq_Max);
	fprintf(parameters, "%8.9g %12.9g %12.9g", params->dist, params->M, params->eta);
	fprintf(parameters, "%12.9g %12.9g %12.9g ", params->bh1.sp, params->bh1.ph, params->bh1.th);
	fprintf(parameters, "%12.9g %12.9g %12.9g ", params->bh2.sp, params->bh2.ph, params->bh2.th);
	fprintf(parameters, "%12.9g %12.9g %12.9g %12.9g ", params->dec, params->pol, params->phi, dt);
	int j;
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%13.9g ", match[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%3d ", good[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(parameters, "%13.9g ", dev[j]);
	}
	fprintf(parameters, "\n");
	fclose(parameters);
}

void write_Data_to_Plotx(ulong i, double match[], double dev[], int good[], int length, const binary_system * const params) {
	FILE * to_Plot;
	if (!(to_Plot = fopen("xto_Plot.data", "at"))) {
		fprintf(stderr, "Can not open the file named \"to_Plot.data\" for append!\n");
		exit(-1);
	}
	fprintf(to_Plot, "%8lu ", i);
	fprintf(to_Plot, "%13.9lg %13.9lg ", params->bh1.m, params->bh2.m);
	fprintf(to_Plot, "%12.9lg %12.9lg ", params->bh1.sp, params->bh2.sp);
	fprintf(to_Plot, "%12.9lg %12.9lg ", params->bh1.ph, params->bh1.th);
	fprintf(to_Plot, "%12.9lg %12.9lg\n", params->bh2.ph, params->bh2.th);
	int j;
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%13.9g ", match[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%3d ", good[j]);
	}
	for (j = 0; j < length; j++) {
		fprintf(to_Plot, "%13.9g ", dev[j]);
	}
	fprintf(to_Plot, "\n");
	fclose(to_Plot);
}
