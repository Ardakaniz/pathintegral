#include "euler_solving.h"

#include <stdio.h>
#include <stdlib.h>

int path_init(path_t* path) {
	const unsigned int path_size = path->N + 1

	path->xs = malloc(path_size * sizeof(double));

	if (path->xs == NULL) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", path_size*sizeof(double));

		return -1;
	}

	path->xps = malloc(path_size * sizeof(double));

	if (path->xps == NULL) {
		fprintf(stderr, "Failed to allocate %lu bytes of memory\n", path_size*sizeof(double));
		free(path->xs);
		return -1;
	}

	return 0;
}

void path_free(path_t* path) {
	if (path == NULL)
		return;

	free(path->xs);
	free(path->xps);
}

void solve_euler_lagrange(path_t* path) {
	for (unsigned int n = 0; n < path->N; ++n) {
		double* x_n = &path->xs[n];
		double* xp_n = &path->xps[n];
		const double dt = path->dt;
		const double t_n = n * dt;

		const double k1 = -Vp(*x_n,                                         t_n) / M;
		const double k2 = -Vp(*x_n + dt / 2.0 * *xp_n,                      t_n) / M;
		const double k3 = -Vp(*x_n + dt / 2.0 * *xp_n + dt * dt / 4.0 * k1, t_n) / M;
		const double k4 = -Vp(*x_n + dt * *xp_n + dt * dt / 2.0 * k2,       t_n) / M;

		*(x_n + 1)  = *x_n + dt * *xp_n + dt * dt / 6 * (k1 + k2 + k3);
		*(xp_n + 1) = *xp_n + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
}

int shoot_and_try(path_t* path) {
	const unsigned int MAX_ITER_COUNT = 25;

	path->xs[0] = path->x_i;
	path->xs[path->N] = path->x_f;
	double xp0_0 = path->xps[0], xp0_1 = path->xps[0];
	
	double err = 0, prev_err = 0;
	for (unsigned int iter_count = 0; iter_count < MAX_ITER_COUNT; ++iter_count) {
		path->xps[0] = xp0_1;
		rk4(path);

		err = fabs(path->xs[path->N] - path->x_f);
		if (err < 1e-3 * SCALE_FACTOR)
			return 0; // Yeah!!, we found the solution
		else { // We did not find the solution yet
			if (iter_count == 0)
				xp0_1 = (2.0 * path->x_f - path->xs[path->N] - path->x_i) / path->t_f;
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
