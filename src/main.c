#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "functions.h"
#include "euler_solving.h"

double Vp(double x, double t) {
	return (V(x + 0.001 * SCALE_FACTOR, t) - V(x, t)) / (0.001 * SCALE_FACTOR);
}

double complex compute_propagator(double* xs, double* xps, double t_f, unsigned int N) {
	const double dt = t_f / (double)N;

	solve_euler_lagrange(xs, xps, N, dt);
	const double action = compute_action(xs, xps, N, dt);
	const double prefactor = compute_prefactor(xs, N, dt);
	
	const double complex propag = cexp(I * action / HBAR) * csqrt(M / (2.0 * PI * I * HBAR * prefactor));

	return propag;
}

double complex compute_wave_fn(double x, double t, double* xs, double* xps, unsigned int N) {
	// How to set x_inf ? c.f. https://www.physics.mcgill.ca/~hilke/719.pdf p.30

	const double x_inf = 20 * sqrt(2 * HBAR * t / M); // <- '20' = arbitrary value...
	const unsigned int p_inf = 150;
	const double dx = x_inf / (double)p_inf;

	double complex partial_wave_fn = 0.0;

	xs[0] = -x_inf + x;
	xs[N] = x;
	partial_wave_fn += compute_propagator(xs, xps, t, N) * wave_fn(-x_inf + x);

	xs[0] = x_inf;
	xs[N] = x;
	partial_wave_fn += compute_propagator(xs, xps, t, N) * wave_fn(x_inf + x);

	partial_wave_fn /= 2.0;

	for (int p = -(int)p_inf + 1; p <= (int)p_inf - 1; ++p) {
		const double x0 = p * dx + x;
		xs[0] = x0;
		xs[N] = x;
		partial_wave_fn += compute_propagator(xs, xps, t, N) * wave_fn(x0);
	}

	return partial_wave_fn * dx;
}

int main(void) {
	const unsigned int N = 100;
	double* xs = malloc((N + 1) * sizeof(double));
	if (xs == NULL) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", (N+1)*sizeof(double));
		
		return EXIT_FAILURE;
	}

	double* xps = malloc((N + 1) * sizeof(double));
	if (xps == NULL) {
		free(xs);

		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", (N+1)*sizeof(double));

		return EXIT_FAILURE;
	}

	FILE* file = fopen("out.dat", "w");

	if (file == NULL) {
		free(xps);
		free(xs);

		fprintf(stderr, "Failed to write to file 'out.dat'\n");

		return EXIT_FAILURE;
	}

	clock_t begin = clock();

	const double x_inf = 0.6 * SCALE_FACTOR;
	const int p_inf = 100;
	const double dx = x_inf / (double)p_inf;

	const double t_min = 0.001 * TIME_FACTOR, t_max = 0.2 * TIME_FACTOR;
	const unsigned int t_count = 100;
	const double dt = (t_max - t_min) / (double)t_count;

	fprintf(file, "%le %le %le %u ", t_min / TIME_FACTOR, t_max / TIME_FACTOR, x_inf / SCALE_FACTOR, (2 * p_inf + 1));

	for (unsigned int i = 0; i < t_count; ++i) {
		if (i % (t_count / 10) == 0) {
			printf("%.0f%%...", 100.0 * (double)i / (double)t_count);
			fflush(stdout);
		}

		const double t = i * dt + t_min;
		
		for (int p = -p_inf; p <= p_inf; ++p) {
			const double x = p * dx;

			const double complex xt_wave_fn = compute_wave_fn(x, t, xs, xps, N);
			const double sq_norm = creal(xt_wave_fn * conj(xt_wave_fn)) * SCALE_FACTOR;

			fprintf(file, "%le ", sq_norm);
		}
	}

	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Total time exceeded: %f\n", time_spent);

	fclose(file);

	free(xps);
	free(xs);
	return EXIT_SUCCESS;
}

/*
int main(void) {
	parameters_t params;
	params->N = 50;

	if (init_parameters(&params) != 0)
		return EXIT_FAILURE;

	// Simulation window:
	double x_min = 0.0, x_max = 0.0;
	double t_min = 0.0, t_max = 0.0;


	return EXIT_SUCCESS
}
*/
