#include "common.h"
#include <iostream>

using namespace std;

void collision()
{
	int z, y, x;
	int zp, yp, xp;
	int k;
	int pindex, index;
	
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
				for(x = xl; x <= xu; x++)
				{
					zp = z - zlg; yp = y - ylg; xp = x - xlg;
					index = zp*NXYg + yp*NXg + xp;
					if(flag[index] == 1)
					{
						pindex = k*NXYZg + index;
						fp[pindex] = f[pindex] - omega*(f[pindex] - feq(rho[index], u[index], v[index], w[index], k)) + Wi[k]*dt*3.0*e[k][0]*Fx;
					}

				}
}
