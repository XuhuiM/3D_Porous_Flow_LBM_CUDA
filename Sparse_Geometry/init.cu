#include "lb.h"
#include "common.h"

void LB_init()
{
	int k;
	double NX;
	double ut, vt, wt, pt;
	double G;

	NX = 90;
	Lx = 1.0;
	Ly = Lx;
	Lz = Lx;
	gama = 1.0;
	rgama = 1.0/gama;
	dx = Lx/NX;
	dt = dx;
	ci = dx/dt;
	rcc = 3.0/ci/ci;
	tau_f = 1.0;
	wf = 1.0/tau_f;
	nu = (tau_f - 0.5)*dt/3;
	U0 = Re*nu/Ly;
//	G = U0*nu*PI*PI*PI/4;
	G = 1.0e-5;
	Fx = G*rgama;
	Fy = 0.0;
	Fz = 0.0;


	sf[0] = 0.0;
	sf[1] = wf; sf[2] = wf;
	sf[3] = sf[5] = sf[7] = wf;
	sf[9] = sf[11] = sf[13] = sf[14] = sf[15] = wf;
	sf[4] = sf[6] = sf[8] = (16.0*tau_f - 8.0)/(8.0*tau_f - 1.0);
	sf[10] = sf[12] = wf;
	sf[16] = sf[17] = sf[18] = sf[4];
	 
	
	printf("tau_f = %f, Ma = %f, G = %e\n", tau_f, sqrt(3.0*rgama)*U0/ci, G);

	ut = 0.0; 
	vt = wt = 0.0;
	pt = 1.0;
	for(k = 0; k < N; k++)
	{
		f[k+0*size]  = FEQ_0(pt,ut,vt,wt,rgama);
		f[k+1*size]  = FEQ_1(pt,ut,vt,wt,rgama);
		f[k+2*size]  = FEQ_2(pt,ut,vt,wt,rgama);
		f[k+3*size]  = FEQ_3(pt,ut,vt,wt,rgama);
		f[k+4*size]  = FEQ_4(pt,ut,vt,wt,rgama);
		f[k+5*size]  = FEQ_5(pt,ut,vt,wt,rgama);
		f[k+6*size]  = FEQ_6(pt,ut,vt,wt,rgama);
		f[k+7*size]  = FEQ_7(pt,ut,vt,wt,rgama);
		f[k+8*size]  = FEQ_8(pt,ut,vt,wt,rgama);
		f[k+9*size]  = FEQ_9(pt,ut,vt,wt,rgama);

		f[k+10*size] = FEQ_10(pt,ut,vt,wt,rgama);
		f[k+11*size] = FEQ_11(pt,ut,vt,wt,rgama);
		f[k+12*size] = FEQ_12(pt,ut,vt,wt,rgama);
		f[k+13*size] = FEQ_13(pt,ut,vt,wt,rgama);
		f[k+14*size] = FEQ_14(pt,ut,vt,wt,rgama);
		f[k+15*size] = FEQ_15(pt,ut,vt,wt,rgama);
		f[k+16*size] = FEQ_16(pt,ut,vt,wt,rgama);
		f[k+17*size] = FEQ_17(pt,ut,vt,wt,rgama);
		f[k+18*size] = FEQ_18(pt,ut,vt,wt,rgama);

	}	
}
