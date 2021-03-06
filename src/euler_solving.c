#include "euler_solving.h"

#include <stdio.h>
#include <stdlib.h>

int path_init(path_t* path) {
	const unsigned int path_size = path->N + 1;

	path->xs = malloc(path_size * sizeof(double));

	if (!path->xs) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", path_size*sizeof(double));

		return -1;
	}

	path->xps = malloc(path_size * sizeof(double));

	if (!path->xps) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", path_size*sizeof(double));
		free(path->xs);
		return -1;
	}

	return 0;
}

void path_free(path_t* path) {
	if (!path)
		return;

	free(path->xs);
	free(path->xps);
}

double compute_lagrangian(double x, double xp, double t) {
	return M * xp * xp / 2.0 - V(x, t);
}

void solve_euler_lagrange(path_t* path, double* action) {
	for (unsigned int n = 0; n < path->N; ++n) {
		double* xs = path->xs;
		double* xps = path->xps;
		const double dt = path->dt;
		const double t_n = n * dt;

		if (action) {
			if (n == 0)
				*action = compute_lagrangian(xs[n], xps[n], t_n) / 2.0;
			else
				*action += compute_lagrangian(xs[n], xps[n], t_n);
		}

		const double k1 = -Vp(xs[n],                                         t_n) / M;
		const double k2 = -Vp(xs[n] + dt / 2.0 * xs[n],                      t_n) / M;
		const double k3 = -Vp(xs[n] + dt / 2.0 * xs[n] + dt * dt / 4.0 * k1, t_n) / M;
		const double k4 = -Vp(xs[n] + dt * xs[n] + dt * dt / 2.0 * k2,       t_n) / M;

		xs[n + 1]  = xs[n]  + dt * xps[n] + dt * dt / 6 * (k1 + k2 + k3);
		xps[n + 1] = xps[n] + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}

	if (action)
		*action = path->dt * (*action + compute_lagrangian(path->xs[path->N], path->xps[path->N], path->t_f) / 2.0);
}

int shoot_and_try(path_t* path, double* action) {
	// See https://www.csun.edu/~lcaretto/me501a/24%20Boundary%20Value%20Problems.pdf
	// Not adapted for steap potentials...
	const unsigned int MAX_ITER_COUNT = 25;

	path->xs[0] = path->x_i;
	path->xs[path->N] = path->x_f;
	double xp0_0 = path->xps[0], xp0_1 = path->xps[0];
	
	double err = 0, prev_err = 0;
	for (unsigned int iter_count = 0; iter_count < MAX_ITER_COUNT; ++iter_count) {
		path->xps[0] = xp0_1;
		
		solve_euler_lagrange(path, action);

		err = fabs(path->xs[path->N] - path->x_f);
		if (err < 1e-3 * SCALE_FACTOR)
			return 0; // Yeah!!, we found the solution
		else { // We did not find the solution yet
			// Finding next xp0 using linear interpolation
			if (iter_count == 0)
				xp0_1 = xp0_0 + (path->x_f - path->xs[path->N]) / path->t_f;
			else {
				const double tmp = xp0_1;
				xp0_1 -= err * (xp0_1 - xp0_0) / (err - prev_err);
				xp0_0 = tmp;
			}
		}

		prev_err = err;
	}

	return -1;
}
