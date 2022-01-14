#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "functions.h"
#include "euler_solving.h"

struct pos_params_t {
	double x_min, x_max;
	unsigned int N;
};
typedef struct pos_params_t pos_params_t;

struct time_params_t {
	double t_min, t_max;
	unsigned int K;
};
typedef struct time_params_t time_params_t;


/*
	*
*/
double* compute_probabilities(pos_params_t* pos_params, time_params_t* time_params, path_t* path) {
	// First: we compute the additional data we need in order to make the calculations
	// How to set x_i_min/max ? c.f. https://www.physics.mcgill.ca/~hilke/719.pdf p.30
	const double x_i_min = -5 * SCALE_FACTOR;
	const double x_i_max = 5 * SCALE_FACTOR;
	const unsigned int P = 200; 
	
	const double dx_i = (x_i_max - x_i_min) / (double)P;
	const double dx_f = (pos_params->x_max - pos_params->x_min) / pos_params->N;

	double t_min = time_params->t_min;
	double time_dt = (time_params->t_max - t_min) / time_params->K;	

	// pos_params actually holds the values of x_f_min and x_f_max
	// We can calculate the adapted x_i_min and x_i_max requiring minimal iteration count using pos_params and time_params

	double* probabilities           = malloc((time_params->K + 1) * (pos_params->N + 1) * sizeof(double));
	double* last_xpN                = malloc((time_params->K + 1) * sizeof(double));
	double complex* partial_wave_fn = malloc((time_params->K + 1) * sizeof(double complex));

	for (unsigned int n = 0; n <= pos_params->N; ++n) { // x_f
		double first_postfactor;

		for (unsigned int p = 0; p <= P; ++p) { // x_i
			path->x_i = x_i_min            + dx_i * p;
			path->x_f = pos_params->x_min  + dx_f * n;

			double action;
			path->t_f = time_params->t_min;
			for (unsigned int k = 0; path->t_f <= time_params->t_max; ++k, path->t_f += time_dt) {
				path->dt = path->t_f / path->N;
				path->xps[0] = (path->x_f - path->x_i) / path->t_f;

				int result = shoot_and_try(path, &action);

				if (result != 0) { // If we failed to find the trajectory
					//if (path->dt < 1.0e-10) { // FAILED TO FIND, RIP...
						printf("Failed to find trajectory, aborting...\n");
						
						free(probabilities);
						free(last_xpN);
						free(partial_wave_fn);

						return NULL;
					//}

					t_min += (k - 1) * path->dt; // We get back to the previous iteration (or decrease t_f_min if k == 0)
					time_dt *= 0.5;        // We improve the timestep 
					k = 0;
				}

				if (1 /* condition to determine */) {
					const double complex postfactor = cexp(I * action / HBAR) * wave_fn(path->x_i);

					if (p == 0) // We assume here P >= 2
						first_postfactor = postfactor; // we store the postfactor here because we cant compute the derivative of the prefactor (we need p = 1 also)
					else {
						const double complex prefactor = csqrt(-M * (path->xps[path->N] - last_xpN[k]) / (2.0 * PI * I * HBAR * dx_i));
						const double complex factor = prefactor * postfactor;
						
						if (p == 1) {
							const double complex prev_factor = prefactor * first_postfactor;
							partial_wave_fn[k] = prev_factor * 0.5 + factor;
						}
						else if (p == P) {
							partial_wave_fn[k] += factor * 0.5;
							partial_wave_fn[k] *= dx_i;

							// wave function is now complete!
							probabilities[k * (pos_params->N + 1) + n] = creal(partial_wave_fn[k] * conj(partial_wave_fn[k])) * SCALE_FACTOR;
						}
						else
							partial_wave_fn[k] += factor;	
					}
					last_xpN[k] = path->xps[path->N];
				}
			}
		}
	}

	free(last_xpN);
	free(partial_wave_fn);
	return probabilities;
}

void output_probabilities(double* probabilities, pos_params_t* pos_params, time_params_t* time_params, const char* filename) {
	FILE* file = fopen(filename, "w");
	if (!file) {
		fprintf(stderr, "Failed to open file %s!", filename);
		return;
	}

	for (unsigned int i = 0; i < (pos_params->N + 1) * (time_params->K + 1); ++i)
		fprintf(file, "%e ", probabilities[i]);

	fclose(file);
}


int main(void) {
	/********************* Simulation parameters **********************/
	pos_params_t pos_params   = { .x_min = 2.0 * SCALE_FACTOR, .x_max = 8.0 * SCALE_FACTOR, .N = 150 };
	time_params_t time_params = { .t_min = 0.5 * TIME_FACTOR, .t_max = 1.5 * TIME_FACTOR, .K = 150 };
	////////////////////////////////////////////////////////////////////

	path_t path;
	path.N = 50;

	if (path_init(&path) != 0)
		return EXIT_FAILURE;

	clock_t clk_begin = clock();

	double* probabilities = compute_probabilities(&pos_params, &time_params, &path);
	if (probabilities == NULL) {
		path_free(&path);
		return EXIT_FAILURE;
	}

	clock_t clk_end = clock();
	double time_spent = (double)(clk_end - clk_begin) / CLOCKS_PER_SEC;
	printf("Total time exceeded: %fs\n", time_spent);

	output_probabilities(probabilities, &pos_params, &time_params, "out.dat");
	
	path_free(&path);
	return EXIT_SUCCESS;
}
