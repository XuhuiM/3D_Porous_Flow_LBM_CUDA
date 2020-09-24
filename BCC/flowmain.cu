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
	int new_step, goon;
	double err = 1.0;
	clock_t time_begin, time_end;

	int device = 0;
    	cudaSetDevice(device);
    	cudaDeviceProp properties;
    	cudaGetDeviceProperties(&properties, device);
    	printf("Lattice Boltzmann Simulation running on: %s\n", properties.name);

	dim3 threads(BX, 1, 1);
  	dim3 threads_BC(BCX, 1, 1);	
	dim3 grid((NX1+BX-1)/BX, NY1, NZ1);
   	dim3 grid_BC_io((NY1+BCX-1)/BCX, NZ1);
	dim3 grid_BC_fb((NX1+BCX-1)/BCX, NZ1);
	dim3 grid_BC_ub((NX1+BCX-1)/BCX, NY1);
	
	srand((unsigned) time(NULL));

	memalloc();
	geo();
	LB_init();
	datasave();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cudaMemcpy(f_dev,  f, sizeof(double)*Q*size, cudaMemcpyHostToDevice);
	cudaMemcpy(flag_dev,  flag, sizeof(int)*size, cudaMemcpyHostToDevice);
	cudaMemcpy(sf_dev, sf, sizeof(double)*Q, cudaMemcpyHostToDevice);
	//////////////////////////////////////////////////////////////////o///////////////////////////////////////////////////////

	time_begin = clock();
loop:
	printf("Enter the num of steps:");
    	scanf("%d", &new_step);
	nmax += new_step;
	printf("nmax = %d\n", nmax);

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

		Evol_flow<<< grid, threads >>>(rgama, sf_dev, dt, Fx, Fy, Fz, fold, fnew);
		Bc_flow_io<<< grid_BC_io, threads_BC >>>(fnew);
		Bc_flow_fb<<< grid_BC_fb, threads_BC >>>(fnew);
		Bc_flow_ub<<< grid_BC_ub, threads_BC >>>(fnew);
		Bc_flow_BB<<< grid, threads >>>(flag_dev, fnew);
	
		if(n%TP == 0)
		{
		       cudaMemcpy(f, f_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost);	
		       err = error();		
		       printf( "n=%d: err = %e\n", n, err);
		       if(n%10000 == 0) datasave();
		}

	}

	time_end = clock();
	printf( "The computing time is: %f mins \n", (double)( time_end - time_begin ) / CLOCKS_PER_SEC / 60.f );
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(n % 2 == 0)
	{
		cudaMemcpy( f, f_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost );
		printf("this is from f\n");
	}
	else
	{
	    cudaMemcpy( f, F_dev, Q*size*sizeof(double), cudaMemcpyDeviceToHost );
	    printf("this is from F\n");
	}

	datasave();

	printf("goon? yes(1) no(0):");
	scanf("%d", &goon);
	if(goon) goto loop;

	
	cudaFree(f_dev);
	cudaFree(F_dev);
	cudaFree(flag_dev);
	cudaFree(sf_dev);
	////////////////////////////////////////////////////////////////////////////////////////////////_GPU
	return 0;

}
