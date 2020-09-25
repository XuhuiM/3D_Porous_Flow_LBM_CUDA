#include "common.h"

void boundary()
{
	int z, y, x;
	int zp, yp, xp;
	int k;
	int index, indexd, pindex, pindexd;
	double Feq, FNeq;
	// for boundary condition at x = xmin
	if(xl == xmin)
	{
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
			{
				zp = z - zlg; yp = y - ylg; xp = xl - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = zp*NXYg + yp*NXg + xp + 1;
				rho[index] = rho[indexd];
				u[index] = u[indexd]; //U0;
				v[index] = v[indexd];
				w[index] = w[indexd];
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;
				}
				
			}
	}

	// for boundary condition at x = xmax
	if(xu == xmax)
	{
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
			{
				zp = z - zlg; yp = y - ylg; xp = xu - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = zp*NXYg + yp*NXg + xp - 1;
				rho[index] = rho[indexd];
				u[index] = u[indexd];
				v[index] = v[indexd];
				w[index] = w[indexd];
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;

				}
				
			}
	}

	// for boundary condition at y = ymin
	if(yl == ymin)
	{
		for(z = zl; z <= zu; z++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = yl - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = zp*NXYg + (yp+1)*NXg + xp;
				rho[index] = rho[indexd];
				u[index] = 0.0;
				v[index] = 0.0;
				w[index] = 0.0;
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;

				}
				
			}
	}

	// for boundary condition at y = ymax
	if(yu == ymax)
	{
		for(z = zl; z <= zu; z++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = yu - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = zp*NXYg + (yp-1)*NXg + xp;
				rho[index] = rho[indexd];
				u[index] = 0.0;
				v[index] = 0.0;
				w[index] = 0.0;
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;

				}
				
			}
	}

	// for boundary condition at z = zmin
	if(zl == zmin)
	{
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = zl - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = (zp+1)*NXYg + yp*NXg + xp;
				rho[index] = rho[indexd];
				u[index] = 0.0;
				v[index] = 0.0;
				w[index] = 0.0;
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;

				}
				
			}
	}
	// for boundary condition at z = zmax
	if(zu == zmax)
	{
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = zu - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				indexd = (zp-1)*NXYg + yp*NXg + xp;
				rho[index] = rho[indexd];
				u[index] = 0.0;
				v[index] = 0.0;
				w[index] = 0.0;
				for(k = 0; k < Q; k++)
				{
					Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
					pindexd = k*NXYZg + indexd;
					FNeq = f[pindexd] - Feq;
					Feq = feq(rho[index], u[index], v[index], w[index], k);
					pindex = k*NXYZg + index;
					f[pindex] = Feq + FNeq;

				}
				
			}
	}

	//z edge
	if((xl == xmin)&&(yl == ymin))
	{
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = yl - ylg; xp = xl - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = zp*NXYg + (yp+1)*NXg + xp + 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}


	}

	if((xl == xmin)&&(yu == ymax))
	{
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = yu - ylg; xp = xl - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = zp*NXYg + (yp-1)*NXg + xp + 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}

	if((xu == xmax)&&(yl == ymin))
	{
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = yl - ylg; xp = xu - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = zp*NXYg + (yp+1)*NXg + xp - 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((xu == xmax)&&(yu == ymax))
	{
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = yu - ylg; xp = xu - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = zp*NXYg + (yp-1)*NXg + xp - 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}

	//y edge
	if((xl == xmin)&&(zl == zmin))
	{
		for(y = yl; y <= yu; y++)
		{
			zp = zl - zlg; yp = y - ylg; xp = xl - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp+1)*NXYg + yp*NXg + xp + 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((xl == xmin)&&(zu == zmax))
	{
		for(y = yl; y <= yu; y++)
		{
			zp = zu - zlg; yp = y - ylg; xp = xl - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp-1)*NXYg + yp*NXg + xp + 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((xu == xmax)&&(zl == zmin))
	{
		for(y = yl; y <= yu; y++)
		{
			zp = zl - zlg; yp = y - ylg; xp = xu - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp+1)*NXYg + yp*NXg + xp - 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((xu == xmax)&&(zu == zmax))
	{
		for(y = yl; y <= yu; y++)
		{
			zp = zu - zlg; yp = y - ylg; xp = xu - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp-1)*NXYg + yp*NXg + xp - 1;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	//x edge
	if((yl == ymin)&&(zl == zmin))
	{
		for(x = xl; x <= xu; x++)
		{
			zp = zl - zlg; yp = yl - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp+1)*NXYg + (yp+1)*NXg + xp;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((yl == ymin)&&(zu == zmax))
	{
		for(x = xl; x <= xu; x++)
		{
			zp = zu - zlg; yp = yl - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp-1)*NXYg + (yp+1)*NXg + xp;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((yu == ymax)&&(zl == zmin))
	{
		for(x = xl; x <= xu; x++)
		{
			zp = zl - zlg; yp = yu - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp+1)*NXYg + (yp-1)*NXg + xp;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}
	if((yu == ymax)&&(zu == zmax))
	{
		for(x = xl; x <= xu; x++)
		{
			zp = zu - zlg; yp = yu - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			indexd = (zp-1)*NXYg + (yp-1)*NXg + xp;
			rho[index] = rho[indexd];
			u[index] = 0.0;
			v[index] = 0.0;
			w[index] = 0.0;
			for(k = 0; k < Q; k++)
			{
				Feq = feq(rho[indexd], u[indexd], v[indexd], w[indexd], k);
				pindexd = k*NXYZg + indexd;
				FNeq = f[pindexd] - Feq;
				Feq = feq(rho[index], u[index], v[index], w[index], k);
				pindex = k*NXYZg + index;
				f[pindex] = Feq + FNeq;

			}
		}

	}	
}
