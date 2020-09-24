#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define NX  90
#define NY  90
#define NZ  90
#define NX1 (NX+1)
#define NY1 (NY+1)
#define NZ1 (NZ+1)
#define Dim  3
#define Q    19
#define size (NZ1*NY1*NX1)
#define GID(z, y, x) (z*NY1*NX1+y*NX1+x)

int e[Q][Dim] = 
{
	{0, 0, 0}, //0

	{1,  0, 0}, //1
	{-1, 0, 0}, //2
	{0,  1, 0}, //3
	{0, -1, 0}, //4
	{0,  0, 1}, //5
	{0,  0,-1}, //6

	{0,  1,  1},//7
	{0, -1, -1},//8
	{0, -1,  1},//9
	{0,  1, -1},//10
	{-1, 0, -1},//11
	{ 1, 0,  1},//12
	{-1, 0,  1},//13
	{ 1, 0, -1},//14
	{-1, 1,  0},//15
	{ 1,-1,  0},//16
	{-1, -1, 0},//17
	{ 1,  1, 0}//18
};


int re[Q] = 
{
	0,
	2,
	1,
	4,
	3,
	6,
	5,


	8,
	7,
	10,
	9,
	12,
	11,
	14,
	13,
	16,
	15,
	18,
	17
};


int NUMF;
int *flag;
int *sparse_id;
int *global2sparse;
int *sparse2global_x, *sparse2global_y, *sparse2global_z;

void tecplot()
{
	int z, y, x, k;
	int q, numf;
	char filename[15];
	FILE *fp;

	sprintf(filename,"%s%s","geo",".plt");
	
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	fprintf(fp,"Title=\"LBM Driven Flow\"\n");
	fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"flag\"\n");
	fprintf(fp,"ZONE I=%d,J=%d,K=%d,F=POINT\n",NX1,NY1,NZ1);
	for(z = 0; z < NZ1; z++)
		for(y = 0; y < NY1; y++)
			for(x = 0; x < NX1; x++)
			{		
				k = GID(z, y, x);
				fprintf(fp,"%d %d %d %d\n", x, y, z, flag[k]);
			}
	fclose(fp);
	

	if((fp=fopen("node_index","w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	fprintf(fp, "%d\n", NUMF);
	for(q = 0; q < Q; q++)
		for(numf = 0; numf < NUMF; numf++)
		{		
			fprintf(fp,"%d ", sparse_id[numf + q*NUMF]);
		}
	fclose(fp);

	if((fp=fopen("index_node","w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	for(numf = 0; numf < NUMF; numf++)
	{		
		fprintf(fp,"%d %d %d\n", sparse2global_x[numf], sparse2global_y[numf], sparse2global_z[numf]);
	}
	fclose(fp);

}

void geo_fluid()
{
	int z, y, x, k;
	int zp, yp, xp, kp;
	int kindx, kindxp;
	int numf = 0;
	int q;

	for(z = 0; z < NZ1; z++)
		for(y = 0; y < NY1; y++)
			for(x = 0; x < NX1; x++)
			{
				k = GID(z, y, x);
				if(flag[k] == 0)
				{
					NUMF += 1;
					global2sparse[k] = numf;
					numf += 1;
				}
			}

	sparse_id = (int*)calloc(NUMF*Q, sizeof(int));
	sparse2global_x = (int*)calloc(NUMF, sizeof(int));
	sparse2global_y = (int*)calloc(NUMF, sizeof(int));
	sparse2global_z = (int*)calloc(NUMF, sizeof(int));

	for(z = 0; z < NZ1;z++)
		for(y = 0; y < NY1; y++)
			for(x = 0; x < NX1; x++)
			{
				k = GID(z, y, x);
				if(flag[k] == 0)
				{
					numf = global2sparse[k];
					sparse2global_x[numf] = x;
					sparse2global_y[numf] = y;
					sparse2global_z[numf] = z;
				}
			}

	for(z = 0; z < NZ1; z++)
		for(y = 0; y < NY1; y++)
			for(x = 0; x < NX1; x++)
			{
				k = GID(z, y, x);
				if(flag[k] == 0)
				{
					kindx = global2sparse[k];
					for(q = 0; q < Q; q++)
					{
						zp = (z + e[q][2] + NZ1)%NZ1;
						yp = (y + e[q][1] + NY1)%NY1; 
						xp = (x + e[q][0] + NX1)%NX1;
						kp = GID(zp, yp, xp);
						if(flag[kp] == 0)
						{
							kindxp = global2sparse[kp];
							sparse_id[kindx + q*NUMF] = kindxp + q*NUMF;
						}
						else
						{
							sparse_id[kindx + q*NUMF] = kindx + re[q]*NUMF;
						}
					}
				}
			}

}


int main()
{
	int z, y, x, k;
	FILE *fp;

	flag = (int*)calloc(size, sizeof(int));
	global2sparse = (int*)calloc(size, sizeof(int));

	if((fp = fopen("fg", "r")) == NULL) {printf("File Open Error\n");exit(1);}
	for(z = 0; z < NZ1; z++)
		for(y = 0; y < NY1; y++)
			for(x = 0; x < NX1; x++)
			{
				k = GID(z, y, x);
				fscanf(fp, "%d ", &flag[k]);
			}
	fclose(fp);


	geo_fluid();

	tecplot();

	free(flag);
	free(global2sparse);
	free(sparse_id);
	return 0;

}
