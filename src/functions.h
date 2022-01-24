#pragma once

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "constants.h"
#include <complex.h>
#include <math.h>

#define SCALE_FACTOR 1e-9
#define TIME_FACTOR 1e-15
#define M 9.1e-31

#define W0 1.0e13

static inline double V(double x, double t) {
	return M * W0 * W0 * x * x / 2.0;
}

static inline double complex wave_fn(double x) {
	return sqrt(sqrt(M * W0 / (PI * HBAR))) * sqrt(2.0 * M * W0 / HBAR) * x * exp(-0.5 * M * W0 / HBAR * x * x);
}

#endif // FUNCTIONS_H