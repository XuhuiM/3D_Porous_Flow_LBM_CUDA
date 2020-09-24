#include "lb.h"
#include "common.h"
#include<stdlib.h>
#include<stdio.h>

void geo()
{
	int k, q;
	FILE *fp;

	if((fp=fopen("node_index", "r"))==NULL) printf("FILE OPEN ERROR!\n");
	fscanf(fp, "%*[^\n]%*c");
	for(q = 0; q < Q; q++)
		for(k = 0; k < N; k++)
		{
			fscanf(fp, "%d ", &node_index[k+q*size]);
		}
	fclose(fp);

	if((fp=fopen("node_inf", "w"))==NULL) printf("FILE OPEN ERROR!\n");
	for(q = 0; q < Q; q++)
		for(k = 0; k < N; k++)
		{
			fprintf(fp, "%d ", node_index[k+q*size]);
		}
	fclose(fp);
}

