extern "C"
{
	void c_cavityflow(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny);
}

#include <stdlib.h>

double **alloc_2d_double(int rows, int cols)
{
	double *data = (double *)malloc(rows*cols*sizeof(int*));
	double **array = (double **)malloc(rows*sizeof(int*));

	for (int i=0; i<rows; i++)
	{
		array[i] = &(data[cols*i]);
	}

	return array;
}

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
	
	c_cavityflow(nt, nit, u, v, dt, dx, dy, p, rho, nu, nx, ny);

  	return 0;
}
