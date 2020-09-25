////input the physical and LBM parameters/////////////////
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "common.h"
#include "mpi.h"

using namespace std;

void parameters()
{
	ifstream infile("input");
    	string para;
	int N;
    	int intParam[11];
    	double dblParam[1];
    	int MPI_ERR;

	if(proc == master)
	{
		cout << "Reading input parameters..." << endl;
		infile >> para >> MaxStep;
		infile >> para >> tStat;   
		infile >> para >> xmin;
        infile >> para >> xmax;
        infile >> para >> ymin;
        infile >> para >> ymax;
        infile >> para >> zmin;
        infile >> para >> zmax;
        infile >> para >> mpi_xdim;
        infile >> para >> mpi_ydim;
        infile >> para >> mpi_zdim;
		infile >> para >> uConv;
        infile.close();
	
		intParam[0]  = MaxStep;
		intParam[1]  = tStat;
		intParam[2]  = xmin;
		intParam[3]  = xmax;
		intParam[4]  = ymin;
		intParam[5]  = ymax;
		intParam[6]  = zmin;
		intParam[7]  = zmax;
		intParam[8]  = mpi_xdim;
		intParam[9]  = mpi_ydim;
		intParam[10] = mpi_zdim;
	
		dblParam[0] = uConv;
	}

	//broadcast input values from master to all other processes
	MPI_Bcast(intParam, 11, MPI_INT, master, MPI_COMM_WORLD);
	MPI_Bcast(dblParam, 1,  MPI_DOUBLE, master, MPI_COMM_WORLD);

	//assign parameters to local-named variables in slave-nodes
	MaxStep  = intParam[0]; 
    	tStat    = intParam[1]; 
    	xmin     = intParam[2]; 
    	xmax     = intParam[3]; 
    	ymin     = intParam[4]; 
    	ymax     = intParam[5]; 
    	zmin     = intParam[6]; 
    	zmax     = intParam[7]; 
    	mpi_xdim = intParam[8];
    	mpi_ydim = intParam[9];
    	mpi_zdim = intParam[10];

    	uConv = dblParam[0];

	Fx = 1.0e-3;
	U0 = 0.01;
    	L = 1.0;
	R = 0.3;
	Lx = L;
	Ly = L;
	Lz = L;
	N = xmax - xmin;
    	dx = L/N;
    	dt = dx;
	tau = 1.0;
	nu = (tau - 0.5)*dt/3.0;
    	Re = U0*L/nu;
    	omega = 1.0/tau;
	if (proc == master) 
        cout << "Re = " << showpoint << Re << " tau = " << tau << " MaxStep = " << MaxStep << endl;		
}
