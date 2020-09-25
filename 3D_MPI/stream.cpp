#include "common.h"

void stream()
{
	int z, y, x;
	int zp, yp, xp;
	int k;
	int zd, yd, xd;
	int index, indexd, indexr;
	int pindex, pindexd, pindexr;

	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
				for(x = xl; x <= xu; x++)
				{
					zp = z - zlg; yp = y - ylg; xp = x - xlg;
					zd = zp - e[k][2]; yd = yp - e[k][1]; xd = xp - e[k][0];
					index   = zp*NXYg + yp*NXg + xp;
					indexd  = zd*NXYg + yd*NXg + xd;
					pindex  = k*NXYZg + index;
					pindexd = k*NXYZg + indexd;
					if((flag[index] == 1)&&(flag[indexd] == 1)) f[pindex] = fp[pindexd];
					else if((flag[index] == 1)&&(flag[indexd] == 0))
					{
						indexr   = zp*NXYg + yp*NXg + xp;
						pindexr  = re[k]*NXYZg + indexr;
						f[pindex] = fp[pindexr];
					}
					else ;	
					
				}
}
