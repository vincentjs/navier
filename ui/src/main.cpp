#include "cfd.h"
#include "cavityFlow.h"

int main(int argc, char* argv[])
{

	const int nx = 41, ny = 41;
	const int nt = 1000, nit = 100;
	int c = 1;

	double dx = 2.0 / (nx-1);
	double dy = 2.0 / (ny-1);

	double rho = 1;
	double nu = 0.1;
	double dt = 0.001;

	double u[nx*ny] = {0};
	double v[nx*ny] = {0};
	double p[nx*ny] = {0};

	CavityFlow solver = CavityFlow();
	
	solver.cavityFlow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);
	
	return 0;
}
