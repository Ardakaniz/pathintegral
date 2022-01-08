#include "functions.h"

inline double Vp(double x, double t) {
	return (V(x + 0.001 * SCALE_FACTOR, t) - V(x, t)) / (0.001 * SCALE_FACTOR);
}