//// computation of the macro variables /////////////////////
#include "common.h"

void relaxstats()
{
	int z, y, x;
	int zp, yp, xp;
	int k;
	int index, pindex;
	double rhot, ut, vt, wt;

	for(z = zl; z <= zu; z++)
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				if(flag[index] == 1)
				{
					rhot = ut = vt = wt = 0.0;
					for(k = 0; k < Q; k++)
					{
						pindex = k*NXYZg + index;
						rhot += f[pindex];
						ut   += e[k][0]*f[pindex];
						vt   += e[k][1]*f[pindex];
						wt   += e[k][2]*f[pindex];
					}
					rho[index] = rhot;
					u[index]   = ut + 0.5*dt*Fx;
					v[index]   = vt;
					w[index]   = wt;
				}
				else rho[index] = u[index] = v[index] = w[index] = 0.0;				
			}
}
