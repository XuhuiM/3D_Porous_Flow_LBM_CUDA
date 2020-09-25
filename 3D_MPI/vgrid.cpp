///////split the whole domain into different subdomains///////////////
#include "mpi.h"
#include "common.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void vgrid()
{
	//local variables
	int complete, direction, partial, shift;
	int px, xeast, xwest, py, ynorth, ysouth, pz, zbottom, ztop;
	int dims[mpi_dim], mpi_coords[mpi_dim];
    	int periodic[mpi_dim];
	int reorder;

	dims[0] = mpi_xdim;
    	dims[1] = mpi_ydim;
    	dims[2] = mpi_zdim;
    	periodic[0] = 1;
    	periodic[1] = 1;
    	periodic[2] = 1;
    	reorder = 1;

	MPI_Cart_create(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, &MPI_COMM_VGRID);
	MPI_Comm_rank(MPI_COMM_VGRID, &vproc);
	MPI_Cart_coords(MPI_COMM_VGRID, vproc, mpi_dim, mpi_coords);

	px = mpi_coords[0];
	py = mpi_coords[1];
	pz = mpi_coords[2];

	//partitioning in the x direction
	complete = (xmax - xmin)/dims[0];
	partial  = (xmax - xmin) - complete*dims[0];
	if((px+1) <= partial)
	{
		xl = xmin + (complete+1)*px;
		xu = xmin + (complete+1)*(px+1) - 1;
	}
	else
	{
		xl = xmin + complete*px + partial;
		xu = xmin + complete*(px+1) + partial - 1;
	}
	if((px+1)%dims[0] == 0) xu = xu + 1;

	//partitioning in the y direction
	complete = (ymax - ymin)/dims[1];
	partial  = (ymax - ymin) - complete*dims[1];
	if((py+1) <= partial)
	{
		yl = ymin + (complete+1)*py;
		yu = ymin + (complete+1)*(py+1) - 1;
	}
	else
	{
		yl = ymin + complete*py + partial;
		yu = ymin + complete*(py+1) + partial - 1;
	}
	if((py+1)%dims[1] == 0) yu = yu + 1;

	//partitioning in the z direction
	complete = (zmax - zmin)/dims[2];
	partial  = (zmax - zmin) - complete*dims[2];
	if((pz+1) <= partial)
	{
		zl = zmin + (complete+1)*pz;
		zu = zmin + (complete+1)*(pz+1) - 1;
	}
	else
	{
		zl = zmin + complete*pz + partial;
		zu = zmin + complete*(pz+1) + partial - 1;
	}
	if((pz+1)%dims[2] == 0) zu = zu + 1;

        //ghost layers	
	xlg = xl - 1;
	ylg = yl - 1;
	zlg = zl - 1;
	xug = xu + 1;
	yug = yu + 1;
	zug = zu + 1;

//	cout << "xl = " << xl << " xu = " << xu << endl;
//	cout << "yl = " << yl << " yu = " << yu << endl;
//	cout << "zl = " << zl << " zu = " << zu << endl;
	
	//determine neighbours of this processor
	direction = 0; // x direction
	shift = 1;
	MPI_Cart_shift(MPI_COMM_VGRID, direction, shift, &west, &east);
	direction = 1; //y direction
//	shift = 1;
	MPI_Cart_shift(MPI_COMM_VGRID, direction, shift, &south, &north);
	direction = 2; //z direction
//	shift = 1;
	MPI_Cart_shift(MPI_COMM_VGRID, direction, shift, &bottom, &top);

	//coordinates of neighbours of this processor in the x/y/z direction
	MPI_Cart_coords(MPI_COMM_VGRID, east, mpi_dim, mpi_coords);
	xeast = mpi_coords[0];
	MPI_Cart_coords(MPI_COMM_VGRID, west, mpi_dim, mpi_coords);
	xwest = mpi_coords[0];
	MPI_Cart_coords(MPI_COMM_VGRID, south, mpi_dim, mpi_coords);
	ysouth = mpi_coords[1];
	MPI_Cart_coords(MPI_COMM_VGRID, north, mpi_dim, mpi_coords);
	ynorth = mpi_coords[1];
	MPI_Cart_coords(MPI_COMM_VGRID, bottom, mpi_dim, mpi_coords);
	zbottom = mpi_coords[2];
	MPI_Cart_coords(MPI_COMM_VGRID, top, mpi_dim, mpi_coords);
	ztop = mpi_coords[2];

	//coordinates of the diagonal neighbours of this processor
	//  Neighbors along x edges
	mpi_coords[0] = px;
	mpi_coords[1] = ynorth;
	mpi_coords[2] = ztop;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &tn);	
	mpi_coords[0] = px;
	mpi_coords[1] = ysouth;
	mpi_coords[2] = ztop;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &ts);	
	mpi_coords[0] = px;
	mpi_coords[1] = ynorth;
	mpi_coords[2] = zbottom;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &bn);	
	mpi_coords[0] = px;
	mpi_coords[1] = ysouth;
	mpi_coords[2] = zbottom;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &bs);
	
	//  Neighbors along y edges
	mpi_coords[0] = xeast;
	mpi_coords[1] = py;
	mpi_coords[2] = ztop;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &te);
	mpi_coords[0] = xwest;
	mpi_coords[1] = py;
   	mpi_coords[2] = ztop;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &tw);
	mpi_coords[0] = xwest;
	mpi_coords[1] = py;
	mpi_coords[2] = zbottom;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &bw);
	mpi_coords[0] = xeast;
	mpi_coords[1] = py;
	mpi_coords[2] = zbottom;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &be);
	
	//  Neighbors along z edges
	mpi_coords[0] = xeast;
	mpi_coords[1] = ynorth;
	mpi_coords[2] = pz;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &ne);
	mpi_coords[0] = xwest;
	mpi_coords[1] = ysouth;
	mpi_coords[2] = pz;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &sw);
	mpi_coords[0] = xwest;
	mpi_coords[1] = ynorth;
	mpi_coords[2] = pz;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &nw);
	mpi_coords[0] = xeast;
	mpi_coords[1] = ysouth;
	mpi_coords[2] = pz;
	MPI_Cart_rank(MPI_COMM_VGRID, mpi_coords, &se);	

//	 cout << " vgrid is done! " << endl;
}
