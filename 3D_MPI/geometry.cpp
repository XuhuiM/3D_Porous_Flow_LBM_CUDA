////// create or read the porous structure /////////////////////
#include "common.h"

void geo()
{
	int z, y, x, k;
	int zp, yp, xp;
	int index;
	int NR = 9;
	double xc, yc, zc;
	double lenx, leny, lenz, l;
	double cc[NR][ndim];

	k = 0;
	for(zc = 0.0; zc <= Lz; zc+=Lz)
		for(yc = 0.0; yc <= Ly; yc+=Ly)
			for(xc = 0.0; xc <= Lx; xc+=Lx)
			{
				cc[k][0] = xc; cc[k][1] = yc; cc[k][2] = zc;
				k++;
			}
	cc[k][0] = Lx/2; cc[k][1] = Ly/2; cc[k][2] = Lz/2;

    for(z = zl; z <= zu; z++)
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				flag[index] = 1;
			}
    for(z = zl; z <= zu; z++)
	    for(y = yl; y <= yu; y++)
		    for(x = xl; x <= xu; x++)
		    {
			    zp = z - zlg; yp = y - ylg; xp = x - xlg;
			    index = zp*NXYg + yp*NXg + xp;
			    lenz = (z - 1)*dx - cc[k][2]; leny = (y - 1)*dx - cc[k][1]; lenx = (x - 1)*dx - cc[k][0];
			    lenz = lenz*lenz; leny = leny*leny; lenx = lenx*lenx;
			    l = lenz + leny + lenx;
			    if(l <= R*R) flag[index] = 0;
		    }
}
