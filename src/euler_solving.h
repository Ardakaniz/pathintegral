#pragma once

#ifndef EULER_SOLVING_H
#define EULER_SOLVING_H

/*
	Structure used to store the data of a path
*/
struct path_t {
	double x_i, x_f, t_f; // Resp. initial position, final position, final time of the given path
	
	double dt;      // Numerical computation constants
	unsigned int N; // Numerical computation constants

	double* xs;  // Positions
	double* xps; // Velocities
};
typedef struct path_t path_t;

/* Allocates path->xs and path->xps with size (N+1) */
int path_init(path_t* path);

/* Free path->xs and path->xps */
void path_free(path_t* path);

// Computes the lagragian at position x, velocity xp and time t
double compute_lagrangian(double x, double xp, double t)

// Computes the derivative of V
inline double Vp(double x, double t);

/*
	4th order runge-kutta method to solve Euler-Lagrange equations with initial conditions
	Initial conditions expected:
		* xs[0]  = initial_position
		* xps[0] = initial_velocity

	Also computes action on-the-fly if action is not NULL
*/
void solve_euler_lagrange(path_t* path, double* action);

/*
	Shoot And Try method to solve Euler-Lagrange equations with boundary conditions
	Boundary conditions expected:
		* xs[0] = initial_position
		* xs[N] = final_position

	Also computes action on-the-fly if action is not NULL
*/
int shoot_and_try(path_t* path, double* action);

#include "euler_solving.inl"

#endif // EULER_SOLVING_H