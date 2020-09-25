////////// Pore-scale LBE solver for NS ////////////////////////////////
/////////  MRT LBE for flow //////////////////////////////////////////////////
///////main function//////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include "common.h"
#include "mpi.h"

using namespace std;

int main(int argc, char **argv)
{

	//initialize MPI enviroment
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	//read and create virtual CPU grid
	parameters();
	vgrid();

	//allocate memory for commom arrays and initialize
	memalloc(1);

	iStep = 0;
	geometry();
	init();
	flagex();
//	datasave(iStep);	
	
	while(iStep < MaxStep)
	{
		iStep++;

		if((proc == master)&&(iStep%tStat == 0)) cout << "iteration step:" << iStep << " " << endl;

		collision();

		postcollision();

		stream();

		relaxstats();
		
//		boundary();

//		if(iStep%tStat == 0) datasave(iStep);

	}
	datasave(iStep);
	memalloc(0);
	MPI_Finalize();
    	return 0;	
}


