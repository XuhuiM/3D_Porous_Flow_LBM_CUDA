#include "lb.h"
#include "common.h"

double error()
{
	int q, k;
	double ut, vt, wt;
	double sum_u_c = 0.0, sum_v_c = 0.0, sum_w_c = 0.0;
	double Err;
	FILE *fp;

	for(k = 0; k < N; k++)
	{
		ut = vt = wt = 0.0;
		for(q = 1; q < Q; q++)
		{
			ut += e[q][0]*f[k+q*size];
			vt += e[q][1]*f[k+q*size];
			wt += e[q][2]*f[k+q*size];

		}	
		ut += 0.5*dt*Fx;
		vt += 0.5*dt*Fy;
		wt += 0.5*dt*Fz;

		sum_u_c += ut;
		sum_v_c += vt;
		sum_w_c += wt;
	}
	

//	Err = fabs(sum_u_c - sum_u_o)/fabs(sum_u_c + 1.0e-9) + fabs(sum_v_c - sum_v_o)/fabs(sum_v_c + 1.0e-9) + fabs(sum_w_c - sum_w_o)/fabs(sum_w_c + 1.0e-9);
	Err = fabs(sum_u_c - sum_u_o)/fabs(sum_u_c + 1.0e-9);
	sum_u_o = sum_u_c;
	sum_v_o = sum_v_c;
	sum_w_o = sum_w_c;

	if((fp=fopen("err.dat", "a"))==NULL) printf("FILE OPEN ERROR!\n");
	fprintf(fp, "%d %e\n", n, Err);
	fclose(fp);

	return Err;

}
