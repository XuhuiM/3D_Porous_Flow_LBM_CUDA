//---LBE for 3D porous flow----------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "cuda_runtime.h"
#include <time.h>
#include "lb.h"
#include "common.h"
#include "memalloc.cu"
#include "geo.cu"
#include "init.cu"
#include "flow.cu"
#include "error.cu"
#include "datasave.cu"
//----------------------------------------------------------------------------------------------------------------------------
int main()
{
	unsigned int new_step, goon;
	cudaEvent_t start, stop;
	float GPU_Time;	
	double err = 1.0;
	clock_t time_begin, time_end;

	int device = 1;
    cudaSetDevice(device);
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, device);
    printf("Lattice Boltzmann Simulation running on: %s\n", properties.name);

	memalloc();

	dim3 threads(BX, 1, 1);
	dim3 grid((N+BX-1)/BX, 1, 1);
	
	geo();
	LB_init();
	datasave();	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cudaMemcpy(f_dev,  f, sizeof(double)*Q*size, cudaMemcpyHostToDevice);
	cudaMemcpy(sf_dev, sf, sizeof(double)*Q, cudaMemcpyHostToDevice);
	cudaMemcpy(node_dev,  node_index, sizeof(int)*Q*size, cudaMemcpyHostToDevice);
	//////////////////////////////////////////////////////////////////o///////////////////////////////////////////////////////

loop:
	printf("Enter the num of steps:");
    scanf("%u", &new_step);
	nmax += new_step;
	printf("nmax = %u\n", nmax);

	time_begin = clock();
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start,0);	
	while((n < nmax)&&(err > 1.0e-6))
	{
	
		n++;
		if(n%2 == 0)
		{
			fold = F_dev;
			fnew = f_dev;
		}
		else
		{
			fold = f_dev;
			fnew = F_dev;
		}

		Evol_flow<<< grid, threads >>>(rgama, sf_dev, dt, Fx, Fy, Fz, N, size, node_dev, fold, fnew);
	
		if(n%TP == 0)
		{
//		       cudaMemcpy(f, f_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost);	
//		       err = error();		
		       printf( "n=%u: err = %e\n", n, err);
//		       datasave();
		}

	}

	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&GPU_Time, start, stop);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	time_end = clock();
	printf("The computing time is: %f s, %f mins\n", (double)( time_end - time_begin ) / CLOCKS_PER_SEC, (double)( time_end - time_begin ) / CLOCKS_PER_SEC / 60.f);
	printf("MUPLS (s) for GPU is %f\n", 1.0*size*n/1000000.0/GPU_Time*1000.0);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(n % 2 == 0)
	{
		cudaMemcpy(f, f_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost);
		printf("this is from f\n");
	}
	else
	{
	    cudaMemcpy(f, F_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost);
	    printf("this is from F\n");
	}

	datasave();

	printf("goon? yes(1) no(0):");
	scanf("%u", &goon);
	if(goon) goto loop;
	
	free(f);
	free(node_index);
	free(sf);
	cudaFree(f_dev);
	cudaFree(F_dev);
	cudaFree(node_dev);
	cudaFree(sf_dev);
	////////////////////////////////////////////////////////////////////////////////////////////////_GPU
	return 0;

}
