#include <stdlib.h>
#include <malloc.h>
#include "common.h"
#include "lb.h"

void memalloc()
{
    //allocate memeory at host
    f    = (double*)calloc(size*Q, sizeof(double));
    nc   = (double*)calloc(NR*Dim, sizeof(double));
    sf   = (double*)calloc(Q, sizeof(double));
    flag = (int*)calloc(size, sizeof(int));

    //allocate memeory at device
    cudaMalloc((void **)&f_dev, size*Q*sizeof(double));
    cudaMalloc((void **)&F_dev, size*Q*sizeof(double));
    cudaMalloc((void **)&sf_dev, Q*sizeof(double));
    cudaMalloc((void **)&flag_dev, size*sizeof(int));
}

