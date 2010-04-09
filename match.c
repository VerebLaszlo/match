/**
 * @file match.c
 * @author László Veréb
 * @date 2010.04.09.
 */

void multi_Malloc(detector * det, size_t len) {
	det->length = len; 
	det->t = fftw_malloc(detector->length * sizeof(double));
	det->s = fftw_malloc(detector->length * sizeof(double));
	det->n = fftw_malloc(detector->length * sizeof(double));
	det->ct = fftw_malloc(detector->length * sizeof(fftw_complex));
	det->cs = fftw_malloc(detector->length * sizeof(fftw_complex));
	det->cn = fftw_malloc(detector->length * sizeof(fftw_complex));
	det->pt = fftw_plan_dft_r2c_1d(detector->length, t, ct, FFTW_ESTIMATE);
	det->ps = fftw_plan_dft_r2c_1d(detector->length, s, cs, FFTW_ESTIMATE);
	det->pn = fftw_plan_dft_r2c_1d(detector->length, n, cn, FFTW_ESTIMATE);
}

void multi_Free(detector det) {
	fftw_free(det.t);
	fftw_free(det.s);
	fftw_free(det.n);
	fftw_free(det.ct);
	fftw_free(det.cs);
	fftw_free(det.cn);
	fftw_destroy_plan(det.pt);
	fftw_destroy_plan(det.ps);
	fftw_destroy_plan(det.pn);
}

void calculate_After(size_t index, detector det[], size_t length) {
	size_t i;
	for (i = 1; i <= length; i++) {
		fftw_execute(det[i]->pt);
		fftw_execute(det[i]->ps);
		fftw_execute(det[i]->pn);
	}
	i++;
	det[i].length = det[length].length;
	det[i].ct = fftw_malloc(det[i].length * sizeof(fftw_complex));
	det[i].cs = fftw_malloc(det[i].length * sizeof(fftw_complex));
	det[i].cn = fftw_malloc(det[i].length * sizeof(fftw_complex));
	det[i].n = fftw_malloc(det[i].length * sizeof(double));
	det[i].ipn = fftw_plan_dft_r2c_1d(det.length, cn, n, FFTW_ESTIMATE);
	if(!index) {
		printf("Correlation: ");
		fflush(stdout);
		calc_Time_Corr(det[i].t, det[i].t, det[i].t, det[i].length);
		printf("h ");
		fflush(stdout);
		calc_Time_Corr(det[i].s, det[i].s, det[i].s, det[i].length);
		printf("s ");
		fflush(stdout);
		calc_Time_Corr(det[i].n, det[i].n, det[i].n, det[i].length);
		printf("n\n");
		fflush(stdout);
	}
}
