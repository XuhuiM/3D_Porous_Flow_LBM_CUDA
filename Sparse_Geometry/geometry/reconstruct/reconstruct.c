#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#define NX  200
#define NY  200
#define NX1 (NX+1)
#define NY1 (NY+1)
#define Dim  2
#define Q    9
#define size (NY1*NX1)
#define GID(y, x)    (y*NX1+x)


double *u, *v, *rho;
int NUMF = 39061;
int *index;
int *index2node_x, *index2node_y;

void tecplot()
{
	int x, y, k;
	int q, numf;
	FILE *fp;

	if((fp=fopen("u","w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	for(y = 0; y < NY1; y++)
	{
		for(x = 0; x < NX1; x++)
		{		
			k = GID(y, x);
			fprintf(fp,"%e ", u[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	if((fp=fopen("v","w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	for(y = 0; y < NY1; y++)
	{
		for(x = 0; x < NX1; x++)
		{		
			k = GID(y, x);
			fprintf(fp,"%e ", v[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	if((fp=fopen("rho","w"))==NULL)
	{
		printf(" File Open Error\n");exit(1);
	}
	for(y = 0; y < NY1; y++)
	{
		for(x = 0; x < NX1; x++)
		{		
			k = GID(y, x);
			fprintf(fp,"%e ", rho[k]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

}

void fluid()
{
	int y, x, k;
	int numf = 0;
	FILE *fp;

	u = (double*)calloc(size, sizeof(double));
	v = (double*)calloc(size, sizeof(double));
	rho = (double*)calloc(size, sizeof(double));
	index2node_x = (int*)calloc(NUMF, sizeof(int));
	index2node_y = (int*)calloc(NUMF, sizeof(int));

	if((fp=fopen("index_node","r"))==NULL) {printf(" File Open Error\n");exit(1);}
	for(numf = 0; numf < NUMF; numf++)
	{		
		fscanf(fp,"%d %d", &index2node_x[numf], &index2node_y[numf]);
	}
	fclose(fp);

	if((fp=fopen("u.dat","r"))==NULL) {printf(" File Open Error\n");exit(1);}
	for(numf = 0; numf < NUMF; numf++)
	{		
		x = index2node_x[numf]; y = index2node_y[numf];
		k = GID(y, x);
		fscanf(fp,"%lf ", &u[k]);
	}
	fclose(fp);

	if((fp=fopen("v.dat","r"))==NULL) {printf(" File Open Error\n");exit(1);}
	for(numf = 0; numf < NUMF; numf++)
	{		
		x = index2node_x[numf]; y = index2node_y[numf];
		k = GID(y, x);
		fscanf(fp,"%lf ", &v[k]);
	}
	fclose(fp);

	if((fp=fopen("rho.dat","r"))==NULL) {printf(" File Open Error\n");exit(1);}
	for(numf = 0; numf < NUMF; numf++)
	{		
		x = index2node_x[numf]; y = index2node_y[numf];
		k = GID(y, x);
		fscanf(fp,"%lf ", &rho[k]);
	}
	fclose(fp);

}


int main()
{
	int y, x, k;

	fluid();

	tecplot();

	return 0;

}
