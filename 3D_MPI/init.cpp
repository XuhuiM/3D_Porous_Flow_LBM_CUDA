#include "common.h"
#include <iostream>

using namespace std;

void init()
{
	int z, y, x;  //global position
	int zp, yp, xp; //processor position
	int k;
	int pindex, index;
	
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
				for(x = xl; x <= xu; x++)
				{
					zp = z - zlg; yp = y - ylg; xp = x - xlg;
					index = zp*NXYg + yp*NXg + xp;
					pindex = k*NXYZg + index;
					if(flag[index] == 1)
					{
						rho[index] = 1.0;
						u[index] = v[index] = w[index] = 0.0;				
						f[pindex] = feq(rho[index], u[index], v[index], w[index], k);
					}
					else f[pindex] = 0.0;
				}
}
