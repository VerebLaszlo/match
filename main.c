#include "match.h"
#include "util.h"
#include "generator.h"
#include "ioHandling.h"

double pi = M_PI;

extern char * program_invocation_short_name;
extern char * program_invocation_name;

typedef struct waves_Tag{
	size_t length;
	double *data;
} waves;

void malloc_Waves(size_t len, waves * wave) {
	wave->length = len;
	wave->data = malloc(wave->length * sizeof(double));
}

int main(int argc, char *argv[]) {
	char mode[30];
	double dt, freq_Min, freq_Max;
	extern binary_system min, max;
	size_t i;
	int j, k, num_Detectors = 0;
	if (argc != 3) {
		printf("The first parameter is the signals maximal length. The second is a file (the file name is max 20 character) with the number of the detectors in the first line, and after that for each detector the next line:\n"
		"detector_name template_file signal_file noise_file");
		fflush(stdout);
	}
	// preinicialization
	size_t max_Length = atoi(argv[1]);
	FILE *file = sfopen_read(argv[2]);												// deallocated line: 124
	fscanf(file, "%d\n", &num_Detectors);
	typedef char type_name[20];
	type_name *t_names = malloc(num_Detectors * sizeof(type_name));					// deallocated line: 153
	type_name *s_names = malloc(num_Detectors * sizeof(type_name));					// deallocated line: 154
	type_name *n_names = malloc(num_Detectors * sizeof(type_name));					// deallocated line: 155
	type_name *d_names = malloc(num_Detectors * sizeof(type_name));					// deallocated line: 174
	for (i = 0; !feof(file); i++) {
		fscanf(file, "%s %s %s %s\n", d_names[i], t_names[i], s_names[i], n_names[i]);
	}
	fclose(file);
	init_generator(mode, &dt, &freq_Min, &freq_Max);
	program_params params = init_program();
	// allocating memory for data and reading it from files
	waves *templates = malloc(num_Detectors * sizeof(waves));						// deallocated line: end
	waves *signals = malloc(num_Detectors * sizeof(waves));							// deallocated line: end
	waves *noises = malloc(num_Detectors * sizeof(waves));							// deallocated line: end
	for (i = 0; i < num_Detectors; i++) {
		malloc_Waves(max_Length, &templates[i]);									// deallocated line: end
		file = sfopen_read(t_names[i]);												// deallocated line: 137
		for (j = 0; !feof(file); j++) {
			fscanf(file, "%*g %lg\n", &templates[i].data[j]);
		}
		templates[i].length = j;
		fclose(file);
		malloc_Waves(max_Length, &signals[i]);										// deallocated line: end
		file = sfopen_read(s_names[i]);												// deallocated line: 144
		for (j = 0; !feof(file); j++) {
			fscanf(file, "%*g %lg\n", &signals[i].data[j]);
		}
		signals[i].length = j;
		fclose(file);
		malloc_Waves(max_Length, &noises[i]);										// deallocated line: end
		file = sfopen_read(n_names[i]);												// deallocated line: 151
		for (j = 0; !feof(file); j++) {
			fscanf(file, "%*g %lg\n", &noises[i].data[j]);
		}
		noises[i].length = j;
		fclose(file);
	}
	free(t_names);
	free(s_names);
	free(n_names);
	// find the shortest length
	detector_Struct *det = malloc(num_Detectors * sizeof(detector_Struct));						// deallocated line: end
	size_t length = fmin(fmin(templates[0].length, signals[0].length),noises[0].length);
	for (i = 1; i < num_Detectors; i++) {
		length = fmin(fmin(templates[i].length, signals[i].length), fmin(noises[i].length, length));
	}
	//	fill the detector structures
	for (i = 0; i < num_Detectors; i++) {
		multi_Malloc(length, &det[i]);												// deallocated line: 221
		det[i].det = con_Det_Str_Enum(d_names[i]);
		for (j = 0; j < length; j++) {
			det[i].t[j] = templates[i].data[j];
			det[i].s[j] = signals[i].data[j];
			det[i].n[j] = noises[i].data[j];
		}
	}
	free(d_names);
	// execute the Fourier algorithm
	for (i = 1; i <= num_Detectors; i++) {
		fftw_execute(det[i].pt);
		fftw_execute(det[i].ps);
		fftw_execute(det[i].pn);
	}
	// calulates the others
	size_t num_Match = num_Detectors;
	calculate_After(0, det, &num_Match);
	// find the frequency band
	double freq_Step, fr = 0.;
	freq_Step = 1. / (dt * length);
	size_t minfr = 0, maxfr = 0;
	while (fr < freq_Min) {
		fr += freq_Step;
		maxfr = ++minfr;
	}
	while (fr < freq_Max) {
		fr += freq_Step;
		maxfr++;
	}
	// calculating the match
	double ts_Sum = 0., tt_Sum = 0., ss_Sum = 0.;
	num_Match++;
	double own_Match[num_Match];
	double match[num_Match ];
	double diff[num_Match];
	int good[num_Match];
	int count[num_Match];
	for (i = 0; i < num_Match - 1; i++) {
		double *norm = psd(det[i].n, det[i].length, dt, blackman);					// deallocated line: 209
		double ts = scalar_freq(det[i].ct, det[i].cs, norm, minfr, maxfr);
		double tt = scalar_freq(det[i].ct, det[i].ct, norm, minfr, maxfr);
		double ss = scalar_freq(det[i].cs, det[i].cs, norm, minfr, maxfr);
		fftw_free(norm);
		if (i < num_Detectors) {
			ts_Sum += ts;
			tt_Sum += tt;
			ss_Sum += ss;
		}
		own_Match[i + 1] = ts / sqrt(tt * ss);
		count[i + 1] = 0;
	}
	own_Match[0] = ts_Sum / sqrt(tt_Sum * ss_Sum);
	count[0] = 0;
	for (i = 0; i < num_Match - 1; i++) {
		multi_Free(det[i]);
	}
	// write the own match to the file
	file = sfopen_write("own.match");												// deallocated line: 231
	for (i = 0; i < num_Match; i++) {
		fprintf(file, "%15.10lg ", own_Match[i]);
	}
	printf("\n");
	fclose(file);
	// Initialize output files
	{
		file = sfopen_write("parameters.data");										// deallocated line: 248
		fprintf(file, "#%7s %13s %13s ", "index", "m1", "m2");
		fprintf(file, "%13s %13s %13s ", "s1x", "s1y", "s1z" );
		fprintf(file, "%13s %13s %13s ", "s2x", "s2y", "s2z" );

		fprintf(file, "%13s %8s %8s ", "inclination", "freq_min", "freq_max" );
		fprintf(file, "%8s %12s %12s ", "distance", "M", "eta" );
		fprintf(file, "%12s %12s %12s ", "s1", "PH1", "TH1" );
		fprintf(file, "%12s %12s %12s ", "s2", "PH2", "TH2" );
		fprintf(file, "%12s %12s %12s %12s", "declination", "polarization", "phi", "dt" );
		
		fprintf(file, "%13s %13s %13s %13s ", "han_match", "liv_match", "hl_match", "cor_match");
		fprintf(file, "%3s %3s %3s %3s ", "h", "l", "hl", "c");
		fprintf(file, "%13s %13s %13s %13s\n", "han_dev", "liv_dev", "hl_dev", "cor_dev");
		fclose(file);
		file = sfopen_write("to_Plot.data");										// deallocated line: 255
		fprintf(file, "#%7s %13s %7s ", "index", "match", "dev");
		fprintf(file, "%13s %13s ", "m1", "m2");
		fprintf(file, "%12s %12s ", "sp1", "sp2");
		fprintf(file, "%12s %12s ", "ph1", "th1");
		fprintf(file, "%12s %12s\n", "ph2", "th2");
		fclose(file);
	}
	srand(time(NULL));
	// starting the loop
	size_t gen, all_m, bad_m, good_m;
	for (i = 0, gen = 0, all_m = 0, bad_m = 0, good_m; ; i++) {
		if (i % 10 == 0) {
			printf("[%d %d]", i, gen);
			fflush(stdout);
		}
		binary_system now = gen_Params();
		//paraméterek száma * számhossz + üreshelyek + módsztring + (fájlnév = 20)
		int string_length = 6 + 12 * 30 + 15 + strlen(mode) + 20;
		char LALparams[string_length];
		sprintf(LALparams, "./lal %30.20lg %30.20lg %30.20lg %30.20lg %30.20lg %30.20lg %30.20lg\
			%30.20lg %30.20lg %3lg %3lg %3lg %30.20lg %s gen.out", now.bh1.m, now.bh2.m,
			now.bh1.sx, now.bh1.sy, now.bh1.sz, now.bh2.sx, now.bh2.sy, now.bh2.sz, now.incl, freq_Min, freq_Max, now.dist, dt, mode);
		system(LALparams);
		errno = 0;
		file = fopen ("gen.out", "r");												// deallocated line: 294
		if (file == NULL) {
			fprintf (stderr, "%s: Couldn't open file for reading %s; %s\n", program_invocation_short_name, "gen.out", strerror (errno));
			fprintf (stderr, "Skipping!\n");
			continue;
		}
		gen++;
		double * fp, * fc;
		fp = malloc(num_Detectors * sizeof(double));											// deallocated line: 295
		fc = malloc(num_Detectors * sizeof(double));											// deallocated line: 296
		for (i = 0; i < num_Detectors; i++) {
			calc_Response_For_Detector(det[i].det, now.dec, now.phi, now.pol, &fp[i], &fc[i]);
		}
		length = 0;
		while (!(feof(file)) && length < max_Length) {
			double p, c;
			fscanf(file, "%*g %lg %lg", &p, &c);
			for (j = 0; j <= num_Detectors; j++) {
				templates[j].data[length] = fp[j] * p + fc[j] * c;
			}
			length++;
		}
		fclose(file);
		free(fp);
		free(fc);
		if (dt * length < params.min_time) {
			continue;
		}
		all_m++;
		length = fmin(fmin(templates[0].length, signals[0].length),noises[0].length);
		for (j = 1; j <= num_Detectors; j++) {
			length = fmin(fmin(templates[j].length, signals[j].length), fmin(noises[j].length, length));
		}
		//	fill the detector structures
		for (j = 0; j < num_Detectors; j++) {
			multi_Malloc(length, &det[j]);												// deallocated line: 353
			//	Már meg van adva, melyik detektor
			for (k = 0; k < length; k++) {
				det[j].t[k] = templates[j].data[k];
				det[j].s[k] = signals[j].data[k];
				det[j].n[k] = noises[j].data[k];
			}
		}
		// execute the Fourier algorithm
		for (j = 1; j <= num_Detectors; j++) {
			fftw_execute(det[j].pt);
			fftw_execute(det[j].ps);
			fftw_execute(det[j].pn);
		}
		// calulates the others
		num_Match = num_Detectors;
		calculate_After(i, det, &num_Match);
		// find the frequency band
		fr = 0.;
		freq_Step = 1. / (dt * length);
		size_t minfr = 0, maxfr = 0;
		while (fr < freq_Min) {
			fr += freq_Step;
			maxfr = ++minfr;
		}
		while (fr < freq_Max) {
			fr += freq_Step;
			maxfr++;
		}
		// calculating the match
		ts_Sum = 0., tt_Sum = 0., ss_Sum = 0.;
		for (j = 0; j < num_Match; j++) {
			double *norm = psd(det[j].n, det[j].length, dt, blackman);					// deallocated line: 343
			double ts = scalar_freq(det[j].ct, det[j].cs, norm, minfr, maxfr);
			double tt = scalar_freq(det[j].ct, det[j].ct, norm, minfr, maxfr);
			double ss = scalar_freq(det[j].cs, det[j].cs, norm, minfr, maxfr);
			fftw_free(norm);
			if (j < num_Detectors) {
				ts_Sum += ts;
				tt_Sum += tt;
				ss_Sum += ss;
			}
			match[j + 1] = ts / sqrt(tt * ss);
		}
		match[0] = ts_Sum / sqrt(tt_Sum * ss_Sum);
		for (j = 0; j < num_Match; j++) {
			multi_Free(det[j]);
		}
		for (j = 0; j < num_Match; j++) {
			double temp1 = own_Match[j] * (1. - params.deviation);
			double temp2 = own_Match[j] * (1. + params.deviation);
			if (temp1 < match[j] && match[j] < temp2) {
				count[j]++;
				good[j] = 1;
				diff[j] = 1. - (match[j] / own_Match[j]);
			} else if (match[j] >= temp2) {
				good[j] = 2;
				diff[j] = 1. - (match[j] / own_Match[j]);
			} else {
				good[j] = 0;
				diff[j] = 1. - (match[j] / own_Match[j]);
			}
		}
		int all_good = 0;
		for (j = 0; j < num_Detectors; j++) {
			all_good += good[j];
		}
		if (all_good == 0) {
			//	FELSZABADíTANI MEMÓRIÁT	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
			continue;
		}
		multi_Malloc(det[0].length * 2, &det[num_Match]);								// deallocated line: 399
		printf("Correlation: ");
		fflush(stdout);
		calc_Time_Corr(det[1].t, det[0].t, det[num_Match].t, det[num_Match].length / 2);
		printf("h ");
		fflush(stdout);
		calc_Time_Corr(det[1].s, det[0].s, det[num_Match].s, det[num_Match].length / 2);
		printf("s ");
		fflush(stdout);
		calc_Time_Corr(det[1].n, det[0].n, det[num_Match].n, det[num_Match].length / 2);
		printf("n\n");
		fflush(stdout);
		fftw_execute(det[num_Match].pt);
		fftw_execute(det[num_Match].ps);
		double *norm = psd(det[num_Match].n, det[num_Match].length, dt, blackman);					// deallocated line: 396
		double ts = scalar_freq(det[num_Match].ct, det[num_Match].cs, norm, minfr, maxfr);
		double tt = scalar_freq(det[num_Match].ct, det[num_Match].ct, norm, minfr, maxfr);
		double ss = scalar_freq(det[num_Match].cs, det[num_Match].cs, norm, minfr, maxfr);
		fftw_free(norm);
		num_Match++;
		match[num_Match] = ts / sqrt(tt * ss);
		multi_Free(det[num_Match]);
		double temp1 = own_Match[num_Match] * (1. - params.deviation);
		double temp2 = own_Match[num_Match] * (1. + params.deviation);
		if (temp1 < match[num_Match] && match[num_Match] < temp2) {
			count[num_Match]++;
			good[num_Match] = 1;
			diff[num_Match] = 1. - (match[num_Match] / own_Match[num_Match]);
		} else if (match[num_Match] >= temp2) {
			good[num_Match] = 2;
			diff[num_Match] = 1. - (match[num_Match] / own_Match[num_Match]);
		} else {
			good[num_Match] = 0;
			diff[num_Match] = 1. - (match[num_Match] / own_Match[num_Match]);
		}
		all_good = 0;
		for (j = 0; j < num_Match; j++) {
			if ( good[j] == 1) {
				all_good = 1;
				break;
			}
		}
		if (all_good == 1) {
			good_m++;
			write_Waveh(i, templates[0].data, length, dt);
			write_Wavel(i, templates[1].data, length, dt);
			write_Gen_Data(i, freq_Min, freq_Max, match, diff, good, &now, dt);
			write_Data_to_Plot(i, match, diff, good, &now);
		} else {
			bad_m++;
			write_Waveh(bad_m, templates[0].data, length, dt);
			write_Wavel(bad_m, templates[1].data, length, dt);
			write_Gen_Datax(bad_m, freq_Min, freq_Max, match, diff, good, &now, dt);
			write_Data_to_Plotx(bad_m, match, diff, good, &now);
		}
		printf("\nindex: %d, gen: %d, all_m: %d, bad_m: %d, good_m: %d\n", i, gen, all_m, bad_m, good_m);
		printf("sum, detecotrs, fft_corr, time_corr: ");
		for (j = 0; j < num_Match; j++) {
			printf("%d ", count[j]);
		}
		printf("\n");
		fflush(stdout);
	}
	// freeing things
	for (i = 0; i < num_Detectors; i++) {
		free(templates[i].data);
		free(signals[i].data);
		free(noises[i].data);
	}
	free(templates);
	free(signals);
	free(noises);
	free(det);
	return 0;
}
