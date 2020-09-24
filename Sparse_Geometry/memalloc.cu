#include <stdlib.h>
#include <malloc.h>
#include "common.h"
#include "lb.h"

void memalloc()
{
	FILE *fp;
	
	if((fp=fopen("node_index", "r"))==NULL) printf("FILE OPEN ERROR!\n");
	fscanf(fp, "%d ", &N);
	fclose(fp);

	printf("Number of fluid node is %d\n", N);
    
	size = N;
	//allocate memeory at host
    f    = (double*)calloc(size*Q, sizeof(double));
    sf   = (double*)calloc(Q, sizeof(double));
    node_index = (int*)calloc(size*Q, sizeof(int));

    //allocate memeory at device
    cudaMalloc((void **)&f_dev, size*Q*sizeof(double));
    cudaMalloc((void **)&F_dev, size*Q*sizeof(double));
    cudaMalloc((void **)&sf_dev, Q*sizeof(double));
    cudaMalloc((void **)&node_dev, size*Q*sizeof(int));
}

