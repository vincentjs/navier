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

/* C Wrappers */
extern "C" CavityFlow* CavityFlow_init()
{
	return new CavityFlow();
}

extern "C" void CavityFlow_cavityFlow(CavityFlow* ptr, int nt, int nit, double u[], int u_sz, double v[], int v_sz, double dt, double dx, double dy, double p[], int p_sz, double rho, double nu, int nx, int ny)
{
	return ptr->cavityFlow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);
}
