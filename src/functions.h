#pragma once

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "constants.h"
#include <complex.h>
#include <math.h>

#define SCALE_FACTOR 1e-9
#define TIME_FACTOR 1e-15
#define M 9.1e-31

inline double V(double x, double t) {
	return 0;
}

inline double complex wave_fn(double x) {
	const double sigma = 0.055 * SCALE_FACTOR;
	const double normalization_factor = 1.0 / sqrt(2 * sigma * sqrt(PI));
	const double dist = 0.5 * SCALE_FACTOR;
	const double p0 = M * 6 * SCALE_FACTOR / TIME_FACTOR;

	//return normalization_factor * exp(-x*x/(2*sigma*sigma));
	return normalization_factor * (exp(-((x + dist) * (x + dist)) / (2.0 * sigma * sigma))*cexp(I * p0 * x / HBAR)
																	+ exp(-((x - dist) * (x - dist)) / (2.0 * sigma * sigma))*cexp(-I * p0 * x / HBAR));
																	//exp(-x*x / (0.5 * sigma * sigma)));//*cexp(I*PI/2));

	

//	return exp(-M * x * x / (2.0 * HBAR)); 
}

#endif // FUNCTIONS_H