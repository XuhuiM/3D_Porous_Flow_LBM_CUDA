#include "common.h"
#include <stdio.h>
#include <stdlib.h>

void datasave(int step)
{
	int z, y, x;
	int zp, yp, xp;
	int index, indexg;
	FILE *filep;

	char filename[50];
	sprintf(filename, "TPL%.8d_%.2d.dat", step, vproc);
	if( (filep = fopen(filename, "w")) == NULL) {printf("FIle Open Error!\n"); exit(-1);}
	fprintf(filep,"Title=\"LBM_MPI FLOW\"\n");
	fprintf(filep,"VARIABLES=\"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"RHO\",\"FLAG\"\n");
	fprintf(filep,"ZONE I=%d,J=%d,K=%d,F=POINT\n",NX,NY,NZ);
	for(z = zl; z <= zu; z++)	
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexg = (z-1)*GNXY + (y-1)*GNX + (x-1);
				rho_o[indexg] = rho[index];
				u_o[indexg]   = u[index];
				v_o[indexg]   = v[index];
				w_o[indexg]   = w[index];
				fprintf(filep, "%d %d %d %e %e %e %e %d\n", x, y, z, u[index], v[index], w[index], rho[index], flag[index]);
			}
	fclose(filep);

}


