#include "lb.h"
#include "common.h"

double error()
{
	int z, y, x, k;
	int q;
	double ut;
	double sum_u_c = 0.0;
	double Err;
	FILE *fp;

	for(z = 1; z <= NZ1; z++)
		for(y = 1; y <= NY1; y++)
			for(x = N16; x < N16+NX1; x++)
			{
				k = GID(z,y,x);
				if(flag[k] == 0)
				{
					ut = 0.0;
					for(q = 1; q < Q; q++)
					{
						ut += e[q][0]*f[k+q*size];
					}	
					ut += 0.5*dt*Fx;
				}
				sum_u_c += ut;
			}
	

	Err = fabs(sum_u_c - sum_u_o)/fabs(sum_u_c);
	sum_u_o = sum_u_c;

	if((fp=fopen("err.dat", "a"))==NULL) printf("FILE OPEN ERROR!\n");
	fprintf(fp, "%d %e\n", n, Err);
	fclose(fp);

	return Err;

}
