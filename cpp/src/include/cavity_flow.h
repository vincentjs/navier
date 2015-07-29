#ifndef CAVITY_FLOW_H
#define CAVITY_FLOW_H

class Cavity_Flow
{
public:
	Cavity_Flow();
	
	void run(int nt, int nit, double u[], double v[], double dt, double dx, double dy, double p[], double rho, double nu, int nx, int ny);

	~Cavity_Flow();
};

#endif

