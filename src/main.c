#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "functions.h"

struct parameters_t {
	double x_i, x_f, t_f;
	double dt;
	unsigned int N;
	

	/*
		General purpose xs and xps arrays to not allocate them numerous times
	*/
	double* xs;
	double* xps;
};
typedef struct parameters_t parameters_t;

int init_parameters(parameters_t* params) {
	params->xs = malloc((params->N + 1) * sizeof(double));

	if (params->xs == NULL) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", (N+1)*sizeof(double));

		return -1;
	}

	params->xps = malloc((params->N + 1) * sizeof(double));

	if (params->xps == NULL) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", (N+1)*sizeof(double));
		free(params->xs);
		return -1;
	}

	return 0;
}

double Vp(double x, double t) {
	return (V(x + 0.001 * SCALE_FACTOR, t) - V(x, t)) / (0.001 * SCALE_FACTOR);
}

// Hopefully temporary
double Vpp(double x, double t) {
	return (Vp(x + 0.001 * SCALE_FACTOR, x) - Vp(x, t)) / (0.001 * SCALE_FACTOR);
}

/*
	* Uses 4th order runge-kutta to solve euler-lagrange equations
	* Initial conditions expected in xs: xs[0] = x_initial ; xs[N] = x_final
*/
void solve_euler_lagrange(parameters_t* params) {
	for (unsigned int n = 0; n < params->N; ++n) {
		double* x_n = &params->xs[n];
		double* xp_n = &params->xps[n];
		const double dt = params->dt;
		const double t_n = n * dt;

		const double k1 = -Vp(*x_n,                                         t_n) / M;
		const double k2 = -Vp(*x_n + dt / 2.0 * *xp_n,                      t_n) / M;
		const double k3 = -Vp(*x_n + dt / 2.0 * *xp_n + dt * dt / 4.0 * k1, t_n) / M;
		const double k4 = -Vp(*x_n + dt * *xp_n + dt * dt / 2.0 * k2,       t_n) / M;

		*(x_n + 1)  = *x_n + dt * *xp_n + dt * dt / 6 * (k1 + k2 + k3);
		*(xp_n + 1) = *xp_n + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}

/*
	L = T - U = m/2 * (dx/dt)² - V(x, t)
*/
double compute_lagrangian(double x_n, double xp_n, double t_n) {
	return M * xp_n * xp_n / 2.0 - V(x_n, t_n);
}

/*
	Computes the action with initial and final positions provided in xs[0] and xs[N] respectively

	Integrates from t_i = 0 to t_n = N * dt

	TODO: No need to integrate when E = cte? c.f. Hamilton-Jacobi
*/
double compute_action(double* xs, double* xps, unsigned int N, double dt) {
	double partial_action = (compute_lagrangian(xs[0], xps[0], 0) + compute_lagrangian(xs[N], xps[N], N * dt)) / 2.0;

	for (unsigned int n = 1; n <= N - 1; ++n)
		partial_action += compute_lagrangian(xs[n], xps[n], n * dt);
	
	return partial_action * dt;
}

/*
	Compute a prefactor for the propagator

	TODO: No need to solve the differential equation, just need to derive the action
*/
double compute_prefactor(double* xs, unsigned int N, double dt) {
	// Initial conditions
	double val   = 0.0;
	double val_p = 1.0;

	// Computation
	for (unsigned int i = 0; i < N; ++i) {
		const double x_n = xs[i];
		const double t_n = i * dt;

		const double k1 = -val / M * Vpp(x_n,                                         t_n);
		const double k2 = -val / M * Vpp(x_n + dt / 2.0 * val_p,                      t_n);
		const double k3 = -val / M * Vpp(x_n + dt / 2.0 * val_p + dt * dt / 4.0 * k1, t_n);
		const double k4 = -val / M * Vpp(x_n + dt * val_p + dt * dt / 2.0 * k2,       t_n);

		val   += dt * val_p + dt * dt / 6 * (k1 + k2 + k3);
		val_p += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}

	return val;
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
