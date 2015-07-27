extern "C"
{
	void c_cavityflow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny);
}

#include <stdlib.h>
#include "cavityFlow.h"

CavityFlow::CavityFlow()
{
}

void CavityFlow::cavityFlow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny)
{
	c_cavityflow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);
}

CavityFlow::~CavityFlow()
{
}
