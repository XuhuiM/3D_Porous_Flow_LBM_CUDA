/////// exchange the distribution functions among different processors ////////////////
#include "mpi.h"
#include "common.h"
#include <iostream>

using namespace std;

void postcollision()
{
	int z, y, x, k;
	int zp, yp, xp;
	int count;
	int index;
	MPI_Request req_snd[18];
    	MPI_Request req_rcv[18];

	//packing
	//x face
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)
			{
				zp = z - zlg; yp = y - ylg; xp = xl - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_west_snd[count] = fp[index];

				zp = z - zlg; yp = y - ylg; xp = xu - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_east_snd[count] = fp[index];

				count++;
			}
	//y face
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = yl - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_south_snd[count] = fp[index];

				zp = z - zlg; yp = yu - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_north_snd[count] = fp[index];

				count++;
			}
	//z face
	count = 0;
	for(k = 0; k < Q; k++)
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)
			{
				zp = zl - zlg; yp = y - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_bot_snd[count] = fp[index];

				zp = zu - zlg; yp = y - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				f_top_snd[count] = fp[index];

				count++;
			}

	//x edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(x = xl; x <= xu; x++)		
		{
			zp = zl - zlg; yp = yl - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_bs_snd[count] = fp[index];

			
			zp = zl - zlg; yp = yu - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_bn_snd[count] = fp[index];

			zp = zu - zlg; yp = yl - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_ts_snd[count] = fp[index];

			zp = zu - zlg; yp = yu - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_tn_snd[count] = fp[index];

			count++;
		}
	//y edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(y = yl; y <= yu; y++)
		{
			zp = zl - zlg; yp = y - ylg; xp = xl - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_bw_snd[count] = fp[index];

			zp = zl - zlg; yp = y - ylg; xp = xu - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_be_snd[count] = fp[index];

			zp = zu - zlg; yp = y - ylg; xp = xl - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_tw_snd[count] = fp[index];

			zp = zu - zlg; yp = y - ylg; xp = xu - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_te_snd[count] = fp[index];

			count++;
		}
	//z edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = yl - ylg; xp = xl - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_sw_snd[count] = fp[index];

			zp = z - zlg; yp = yl - ylg; xp = xu - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_se_snd[count] = fp[index];

			zp = z - zlg; yp = yu - ylg; xp = xl - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_nw_snd[count] = fp[index];

			zp = z - zlg; yp = yu - ylg; xp = xu - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			f_ne_snd[count] = fp[index];

			count++;
		}

	MPI_Irecv(f_west_rcv , xDVsize, MPI_DOUBLE, west , 0, MPI_COMM_VGRID, &req_rcv[0]);
   	MPI_Irecv(f_east_rcv , xDVsize, MPI_DOUBLE, east , 1, MPI_COMM_VGRID, &req_rcv[1]);
    	MPI_Irecv(f_south_rcv, yDVsize, MPI_DOUBLE, south, 2, MPI_COMM_VGRID, &req_rcv[2]);
    	MPI_Irecv(f_north_rcv, yDVsize, MPI_DOUBLE, north, 3, MPI_COMM_VGRID, &req_rcv[3]);
    	MPI_Irecv(f_bot_rcv  , zDVsize, MPI_DOUBLE, bottom, 4, MPI_COMM_VGRID, &req_rcv[4]);
    	MPI_Irecv(f_top_rcv  , zDVsize, MPI_DOUBLE, top  , 5, MPI_COMM_VGRID, &req_rcv[5]);

    	MPI_Irecv(f_bs_rcv   , xedgeDV, MPI_DOUBLE, bs   , 6, MPI_COMM_VGRID, &req_rcv[6]);
    	MPI_Irecv(f_bn_rcv   , xedgeDV, MPI_DOUBLE, bn   , 7, MPI_COMM_VGRID, &req_rcv[7]);
    	MPI_Irecv(f_ts_rcv   , xedgeDV, MPI_DOUBLE, ts   , 8, MPI_COMM_VGRID, &req_rcv[8]);
    	MPI_Irecv(f_tn_rcv   , xedgeDV, MPI_DOUBLE, tn   , 9, MPI_COMM_VGRID, &req_rcv[9]);

    	MPI_Irecv(f_bw_rcv   , yedgeDV, MPI_DOUBLE, bw   , 10, MPI_COMM_VGRID, &req_rcv[10]);
    	MPI_Irecv(f_be_rcv   , yedgeDV, MPI_DOUBLE, be   , 11, MPI_COMM_VGRID, &req_rcv[11]);
    	MPI_Irecv(f_tw_rcv   , yedgeDV, MPI_DOUBLE, tw   , 12, MPI_COMM_VGRID, &req_rcv[12]);
    	MPI_Irecv(f_te_rcv   , yedgeDV, MPI_DOUBLE, te   , 13, MPI_COMM_VGRID, &req_rcv[13]);

    	MPI_Irecv(f_sw_rcv   , zedgeDV, MPI_DOUBLE, sw   , 14, MPI_COMM_VGRID, &req_rcv[14]);
    	MPI_Irecv(f_se_rcv   , zedgeDV, MPI_DOUBLE, se   , 15, MPI_COMM_VGRID, &req_rcv[15]);
    	MPI_Irecv(f_nw_rcv   , zedgeDV, MPI_DOUBLE, nw   , 16, MPI_COMM_VGRID, &req_rcv[16]);
    	MPI_Irecv(f_ne_rcv   , zedgeDV, MPI_DOUBLE, ne   , 17, MPI_COMM_VGRID, &req_rcv[17]);

    	MPI_Isend(f_west_snd , xDVsize, MPI_DOUBLE, west , 1, MPI_COMM_VGRID, &req_snd[0]);
    	MPI_Isend(f_east_snd , xDVsize, MPI_DOUBLE, east , 0, MPI_COMM_VGRID, &req_snd[1]);
   	MPI_Isend(f_south_snd, yDVsize, MPI_DOUBLE, south, 3, MPI_COMM_VGRID, &req_snd[2]);
    	MPI_Isend(f_north_snd, yDVsize, MPI_DOUBLE, north, 2, MPI_COMM_VGRID, &req_snd[3]);
    	MPI_Isend(f_bot_snd  , zDVsize, MPI_DOUBLE, bottom, 5, MPI_COMM_VGRID, &req_snd[4]);
    	MPI_Isend(f_top_snd  , zDVsize, MPI_DOUBLE, top  , 4, MPI_COMM_VGRID, &req_snd[5]);

    	MPI_Isend(f_bs_snd   , xedgeDV, MPI_DOUBLE, bs   , 9, MPI_COMM_VGRID, &req_snd[6]);
    	MPI_Isend(f_bn_snd   , xedgeDV, MPI_DOUBLE, bn   , 8, MPI_COMM_VGRID, &req_snd[7]);
    	MPI_Isend(f_ts_snd   , xedgeDV, MPI_DOUBLE, ts   , 7, MPI_COMM_VGRID, &req_snd[8]);
    	MPI_Isend(f_tn_snd   , xedgeDV, MPI_DOUBLE, tn   , 6, MPI_COMM_VGRID, &req_snd[9]);

    	MPI_Isend(f_bw_snd   , yedgeDV, MPI_DOUBLE, bw   , 13, MPI_COMM_VGRID, &req_snd[10]);
    	MPI_Isend(f_be_snd   , yedgeDV, MPI_DOUBLE, be   , 12, MPI_COMM_VGRID, &req_snd[11]);
    	MPI_Isend(f_tw_snd   , yedgeDV, MPI_DOUBLE, tw   , 11, MPI_COMM_VGRID, &req_snd[12]);
    	MPI_Isend(f_te_snd   , yedgeDV, MPI_DOUBLE, te   , 10, MPI_COMM_VGRID, &req_snd[13]);

    	MPI_Isend(f_sw_snd   , zedgeDV, MPI_DOUBLE, sw   , 17, MPI_COMM_VGRID, &req_snd[14]);
    	MPI_Isend(f_se_snd   , zedgeDV, MPI_DOUBLE, se   , 16, MPI_COMM_VGRID, &req_snd[15]);
    	MPI_Isend(f_nw_snd   , zedgeDV, MPI_DOUBLE, nw   , 15, MPI_COMM_VGRID, &req_snd[16]);
    	MPI_Isend(f_ne_snd   , zedgeDV, MPI_DOUBLE, ne   , 14, MPI_COMM_VGRID, &req_snd[17]);

    	MPI_Waitall(18, req_rcv, MPI_STATUSES_IGNORE);
    	MPI_Waitall(18, req_snd, MPI_STATUSES_IGNORE);

	//unpacking
	//x face	
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(y = yl; y <= yu; y++)			
			{
				
				zp = z - zlg; yp = y - ylg; xp = xug - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_east_rcv[count];

				zp = z - zlg; yp = y - ylg; xp = xlg - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_west_rcv[count];

				count++;
			}
	//y face
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
			for(x = xl; x <= xu; x++)
			{
				zp = z - zlg; yp = ylg - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_south_rcv[count];

				zp = z - zlg; yp = yug - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_north_rcv[count];

				count++;
			}	
	//z face
	count = 0;
	for(k = 0; k < Q; k++)
		for(y = yl; y <= yu; y++)
			for(x = xl; x <= xu; x++)			
			{
				zp = zlg - zlg; yp = y - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_bot_rcv[count];

				zp = zug - zlg; yp = y - ylg; xp = x - xlg;
				index = k*NXYZg + zp*NXYg + yp*NXg + xp;
				fp[index] = f_top_rcv[count];

				count++;
			}
	//x edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(x = xl; x <= xu; x++)
		{
			zp = zlg - zlg; yp = ylg - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_bs_rcv[count];

			zp = zlg - zlg; yp = yug - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_bn_rcv[count];

			zp = zug - zlg; yp = ylg - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_ts_rcv[count];

			zp = zug - zlg; yp = yug - ylg; xp = x - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_tn_rcv[count];

			count++;
		}
	//y edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(y = yl; y <= yu; y++)		
		{
			zp = zlg - zlg; yp = y - ylg; xp = xlg - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_bw_rcv[count];

			zp = zlg - zlg; yp = y - ylg; xp = xug - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_be_rcv[count];

			zp = zug - zlg; yp = y - ylg; xp = xlg - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_tw_rcv[count];

			zp = zug - zlg; yp = y - ylg; xp = xug - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_te_rcv[count];

			count++;
		}
	//z edge
	count = 0;
	for(k = 0; k < Q; k++)
		for(z = zl; z <= zu; z++)
		{
			zp = z - zlg; yp = ylg - ylg; xp = xlg - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_sw_rcv[count];

			zp = z - zlg; yp = ylg - ylg; xp = xug - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_se_rcv[count];

			zp = z - zlg; yp = yug - ylg; xp = xlg - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_nw_rcv[count];

			zp = z - zlg; yp = yug - ylg; xp = xug - xlg;
			index = k*NXYZg + zp*NXYg + yp*NXg + xp;
			fp[index] = f_ne_rcv[count];

			count++;
		}
}
