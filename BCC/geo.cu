#include "lb.h"
#include "common.h"
#include<stdlib.h>
#include<stdio.h>

void geo()
{

	int z, y, x, k, kk;
	int zp, yp, xp, kp;
	int zd, yd, xd, kd;
	int q;
	double zc, yc, xc;
	double lx, ly, lz, l;
	int zi;
	FILE *fp;

	dx = Lx/NX;

	if((fp = fopen("fg", "r"))==NULL) printf("FILE OPEN ERROR!\n");
	for(z=1; z<=NZ1; z++)
		for(y=1; y<=NY1; y++) 
                for(x=N16; x<(NX1+N16); x++)
				{
					k = GID(z, y, x);
					fscanf(fp, "%d ", &flag[k]);
			
				}
	fclose(fp);

	//inlet
	x = N16;
	for(z = 1; z <= NZ1; z++)
		for(y = 1; y <= NY1; y++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(xp < N16) xp = N16+NX;
				else;
				if(yp < 1)   yp = NY1;
				else if(yp > NY1) yp = 1;
				else;
				if(zp < 1)   zp = NZ1;
				else if(zp > NZ1) zp = 1;
				else;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}
	//outlet
	x = N16+NX;
	for(z = 1; z <= NZ1; z++)
		for(y = 1; y <= NY1; y++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(xp > N16+NX) xp = 1;
				else;
				if(yp < 1)   yp = NY1;
				else if(yp > NY1) yp = 1;
				else;
				if(zp < 1)   zp = NZ1;
				else if(zp > NZ1) zp = 1;
				else;
				if(xp > N16+NX) xp = N16;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}

	//front
	y = 1;
	for(z = 1; z <= NZ1; z++)
		for(x = N16+1; x < N16+NX; x++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(yp < 1) yp = NY1;
				else;
				if(zp < 1)   zp = NZ1;
				else if(zp > NZ1) zp = 1;
				else;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}
	//back
	y = NY1;
	for(z = 1; z <= NZ1; z++)
		for(x = N16+1; x < N16+NX; x++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(yp > NY1) yp = 1;
				else;
				if(zp < 1)   zp = NZ1;
				else if(zp > NZ1) zp = 1;
				else;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}
	//bottom
	z = 1;
	for(y = 2; y < NY1; y++)
		for(x = N16+1; x < N16+NX; x++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(zp < 1) zp = NZ1;
				else;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}

	//upper
	z = NZ1;
	for(y = 2; y < NY1; y++)
		for(x = N16+1; x < N16+NX; x++)
			for(q = 1; q < Q; q++)
			{
				xp = x - e[q][0]; yp = y - e[q][1]; zp = z - e[q][2];
				xd = xp; yd = yp; zd = zp;
				if(zp > NZ1) zp = 1;
				else;
				kp = GID(zp, yp, xp);
				if(flag[kp] == 1)
				{
					kd = GID(zd, yd, xd);
					flag[kd] = 1;

				}
			}

	
	zi = NZ2/2;
	if((fp = fopen("flag_zc.dat", "w"))==NULL) printf("FILE OPEN ERROR!\n");
	for(y=1; y<=NY1; y++) 
	{
		for(x=N16; x<(NX1+N16); x++)
		{
			k = GID(zi, y, x);
			fprintf(fp, "%d ", flag[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	zi = 1;
	if((fp = fopen("flag_zi.dat", "w"))==NULL) printf("FILE OPEN ERROR!\n");
	for(y=1; y<=NY1; y++) 
	{
		for(x=N16; x<(NX1+N16); x++)
		{
			k = GID(zi, y, x);
			fprintf(fp, "%d ", flag[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

}

