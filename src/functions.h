#pragma once

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "constants.h"
#include <complex.h>
#include <math.h>

#define SCALE_FACTOR 1e-9
#define TIME_FACTOR 1e-15
#define M 9.1e-31

static inline double V(double x, double t) {
	return 0.0;
}

static inline double complex wave_fn(double x) {
	const double sigma = 0.055 * SCALE_FACTOR;
	const double p0 = M * 6 * SCALE_FACTOR / TIME_FACTOR;
	const double normalization_factor = 1.0 / sqrt(sqrt(PI * sigma * sigma));

	return normalization_factor * cexp(-x * x / (2.0 * sigma * sigma)) * cexp(I * x * p0 / HBAR);
}

#endif // FUNCTIONS_H