extern "C"
{
	void c_cavityflow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny);
}

#include <cavity_flow.h>

Cavity_Flow::Cavity_Flow()
{
}

/*void Cavity_Flow::run(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny)
{
	c_cavityflow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);
	}*/

Cavity_Flow::~Cavity_Flow()
{
}

