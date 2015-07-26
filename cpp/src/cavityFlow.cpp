extern "C"
{
	void c_cavityflow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny);
}

#include <stdlib.h>
#include "cavityFlow.h"

void cavityFlow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny)
{
	/*const int nx = 41, ny = 41;
	const int nt = 1000, nit = 100;
	int c = 1;

	double dx = 2.0 / (nx-1);
	double dy = 2.0 / (ny-1);

	double rho = 1;
	double nu = 0.1;
	double dt = 0.001;

	double f_u[nx*ny] = {0};
	double f_v[nx*ny] = {0};
	double f_p[nx*ny] = {0};*/
	
	c_cavityflow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);
}
