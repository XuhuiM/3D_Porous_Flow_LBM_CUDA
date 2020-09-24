#include "lb.h"
#include "common.h"


// save the results
void datasave()
{
	int q, k;
	double rhot, ut, vt, wt;
	FILE *fp1, *fp2, *fp3, *fp4;

	if((fp1=fopen("u.dat","w")) == NULL) return;
	if((fp2=fopen("v.dat","w")) == NULL) return;
	if((fp3=fopen("w.dat","w")) == NULL) return;
	if((fp4=fopen("rho.dat","w")) == NULL) return;

	for(k = 0; k < N; k++)
	{
		ut = vt = wt = rhot = 0.0;
		for(q = 0; q < Q; q++)
		{
			ut   += e[q][0]*f[k+q*size];
			vt   += e[q][1]*f[k+q*size];
			wt   += e[q][2]*f[k+q*size];
			rhot += f[k+q*size];
		}
		ut += 0.5*dt*Fx;
		vt += 0.5*dt*Fy;
		wt += 0.5*dt*Fz;

		fprintf(fp1, "%e ", ut);
		fprintf(fp2, "%e ", vt);
		fprintf(fp3, "%e ", wt);
		fprintf(fp4, "%e ", rhot);
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
