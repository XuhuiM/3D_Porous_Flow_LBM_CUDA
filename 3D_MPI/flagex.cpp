#include "mpi.h"
#include "common.h"

void flagex()
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
	for(z = zl; z <= zu; z++)
		for(y = yl; y <= yu; y++)
		{
			zp = z - zlg; yp = y - ylg; xp = xl - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_west_snd[count] = flag[index];

			zp = z - zlg; yp = y - ylg; xp = xu - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_east_snd[count] = flag[index];

			count++;
		}
	//y face
	count = 0;
	for(z = zl; z <= zu; z++)
		for(x = xl; x <= xu; x++)
		{
			zp = z - zlg; yp = yl - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_south_snd[count] = flag[index];

			zp = z - zlg; yp = yu - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_north_snd[count] = flag[index];

			count++;
		}
	//z face
	count = 0;
	for(y = yl; y <= yu; y++)
		for(x = xl; x <= xu; x++)
		{
			zp = zl - zlg; yp = y - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_bot_snd[count] = flag[index];

			zp = zu - zlg; yp = y - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag_top_snd[count] = flag[index];

			count++;
		}

	//x edge
	count = 0;
	for(x = xl; x <= xu; x++)		
	{
		zp = zl - zlg; yp = yl - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_bs_snd[count] = flag[index];

			
		zp = zl - zlg; yp = yu - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_bn_snd[count] = flag[index];

		zp = zu - zlg; yp = yl - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_ts_snd[count] = flag[index];

		zp = zu - zlg; yp = yu - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_tn_snd[count] = flag[index];

		count++;
	}
	//y edge
	count = 0;
	for(y = yl; y <= yu; y++)
	{
		zp = zl - zlg; yp = y - ylg; xp = xl - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_bw_snd[count] = flag[index];

		zp = zl - zlg; yp = y - ylg; xp = xu - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_be_snd[count] = flag[index];

		zp = zu - zlg; yp = y - ylg; xp = xl - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_tw_snd[count] = flag[index];

		zp = zu - zlg; yp = y - ylg; xp = xu - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_te_snd[count] = flag[index];

		count++;
	}
	//z edge
	count = 0;
	for(z = zl; z <= zu; z++)
	{
		zp = z - zlg; yp = yl - ylg; xp = xl - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_sw_snd[count] = flag[index];

		zp = z - zlg; yp = yl - ylg; xp = xu - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_se_snd[count] = flag[index];

		zp = z - zlg; yp = yu - ylg; xp = xl - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_nw_snd[count] = flag[index];

		zp = z - zlg; yp = yu - ylg; xp = xu - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag_ne_snd[count] = flag[index];

		count++;
	}

	MPI_Irecv(flag_west_rcv , xsize, MPI_INT, west , 0, MPI_COMM_VGRID, &req_rcv[0]);
	MPI_Irecv(flag_east_rcv , xsize, MPI_INT, east , 1, MPI_COMM_VGRID, &req_rcv[1]);
    	MPI_Irecv(flag_south_rcv, ysize, MPI_INT, south, 2, MPI_COMM_VGRID, &req_rcv[2]);
    	MPI_Irecv(flag_north_rcv, ysize, MPI_INT, north, 3, MPI_COMM_VGRID, &req_rcv[3]);
    	MPI_Irecv(flag_bot_rcv  , zsize, MPI_INT, bottom, 4, MPI_COMM_VGRID, &req_rcv[4]);
    	MPI_Irecv(flag_top_rcv  , zsize, MPI_INT, top  , 5, MPI_COMM_VGRID, &req_rcv[5]);

    	MPI_Irecv(flag_bs_rcv   , xedge, MPI_INT, bs   , 6, MPI_COMM_VGRID, &req_rcv[6]);
    	MPI_Irecv(flag_bn_rcv   , xedge, MPI_INT, bn   , 7, MPI_COMM_VGRID, &req_rcv[7]);
    	MPI_Irecv(flag_ts_rcv   , xedge, MPI_INT, ts   , 8, MPI_COMM_VGRID, &req_rcv[8]);
    	MPI_Irecv(flag_tn_rcv   , xedge, MPI_INT, tn   , 9, MPI_COMM_VGRID, &req_rcv[9]);

    	MPI_Irecv(flag_bw_rcv   , yedge, MPI_INT, bw   , 10, MPI_COMM_VGRID, &req_rcv[10]);
    	MPI_Irecv(flag_be_rcv   , yedge, MPI_INT, be   , 11, MPI_COMM_VGRID, &req_rcv[11]);
    	MPI_Irecv(flag_tw_rcv   , yedge, MPI_INT, tw   , 12, MPI_COMM_VGRID, &req_rcv[12]);
    	MPI_Irecv(flag_te_rcv   , yedge, MPI_INT, te   , 13, MPI_COMM_VGRID, &req_rcv[13]);

    	MPI_Irecv(flag_sw_rcv   , zedge, MPI_INT, sw   , 14, MPI_COMM_VGRID, &req_rcv[14]);
    	MPI_Irecv(flag_se_rcv   , zedge, MPI_INT, se   , 15, MPI_COMM_VGRID, &req_rcv[15]);
    	MPI_Irecv(flag_nw_rcv   , zedge, MPI_INT, nw   , 16, MPI_COMM_VGRID, &req_rcv[16]);
    	MPI_Irecv(flag_ne_rcv   , zedge, MPI_INT, ne   , 17, MPI_COMM_VGRID, &req_rcv[17]);

    	MPI_Isend(flag_west_snd , xsize, MPI_INT, west , 1, MPI_COMM_VGRID, &req_snd[0]);
    	MPI_Isend(flag_east_snd , xsize, MPI_INT, east , 0, MPI_COMM_VGRID, &req_snd[1]);
	MPI_Isend(flag_south_snd, ysize, MPI_INT, south, 3, MPI_COMM_VGRID, &req_snd[2]);
    	MPI_Isend(flag_north_snd, ysize, MPI_INT, north, 2, MPI_COMM_VGRID, &req_snd[3]);
    	MPI_Isend(flag_bot_snd  , zsize, MPI_INT, bottom, 5, MPI_COMM_VGRID, &req_snd[4]);
    	MPI_Isend(flag_top_snd  , zsize, MPI_INT, top  , 4, MPI_COMM_VGRID, &req_snd[5]);

    	MPI_Isend(flag_bs_snd   , xedge, MPI_INT, bs   , 9, MPI_COMM_VGRID, &req_snd[6]);
    	MPI_Isend(flag_bn_snd   , xedge, MPI_INT, bn   , 8, MPI_COMM_VGRID, &req_snd[7]);
    	MPI_Isend(flag_ts_snd   , xedge, MPI_INT, ts   , 7, MPI_COMM_VGRID, &req_snd[8]);
    	MPI_Isend(flag_tn_snd   , xedge, MPI_INT, tn   , 6, MPI_COMM_VGRID, &req_snd[9]);

    	MPI_Isend(flag_bw_snd   , yedge, MPI_INT, bw   , 13, MPI_COMM_VGRID, &req_snd[10]);
    	MPI_Isend(flag_be_snd   , yedge, MPI_INT, be   , 12, MPI_COMM_VGRID, &req_snd[11]);
    	MPI_Isend(flag_tw_snd   , yedge, MPI_INT, tw   , 11, MPI_COMM_VGRID, &req_snd[12]);
    	MPI_Isend(flag_te_snd   , yedge, MPI_INT, te   , 10, MPI_COMM_VGRID, &req_snd[13]);

    	MPI_Isend(flag_sw_snd   , zedge, MPI_INT, sw   , 17, MPI_COMM_VGRID, &req_snd[14]);
    	MPI_Isend(flag_se_snd   , zedge, MPI_INT, se   , 16, MPI_COMM_VGRID, &req_snd[15]);
    	MPI_Isend(flag_nw_snd   , zedge, MPI_INT, nw   , 15, MPI_COMM_VGRID, &req_snd[16]);
    	MPI_Isend(flag_ne_snd   , zedge, MPI_INT, ne   , 14, MPI_COMM_VGRID, &req_snd[17]);

    	MPI_Waitall(18, req_rcv, MPI_STATUSES_IGNORE);
    	MPI_Waitall(18, req_snd, MPI_STATUSES_IGNORE);

	//unpacking
	//x face	
	count = 0;
	for(z = zl; z <= zu; z++)
		for(y = yl; y <= yu; y++)			
		{
				
			zp = z - zlg; yp = y - ylg; xp = xug - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag[index] = flag_east_rcv[count];

			zp = z - zlg; yp = y - ylg; xp = xlg - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag[index] = flag_west_rcv[count];

			count++;
		}
	//y face
	count = 0;
	for(z = zl; z <= zu; z++)
		for(x = xl; x <= xu; x++)
		{
			zp = z - zlg; yp = ylg - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag[index] = flag_south_rcv[count];

			zp = z - zlg; yp = yug - ylg; xp = x - xlg;
			index = zp*NXYg + yp*NXg + xp;
			flag[index] = flag_north_rcv[count];

			count++;
		}	
	//z face
	count = 0;
	for(y = yl; y <= yu; y++)
		for(x = xl; x <= xu; x++)			
			{
				zp = zlg - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				flag[index] = flag_bot_rcv[count];

				zp = zug - zlg; yp = y - ylg; xp = x - xlg;
				index = zp*NXYg + yp*NXg + xp;
				flag[index] = flag_top_rcv[count];

				count++;
			}
	//x edge
	count = 0;
	for(x = xl; x <= xu; x++)
	{
		zp = zlg - zlg; yp = ylg - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_bs_rcv[count];

		zp = zlg - zlg; yp = yug - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_bn_rcv[count];

		zp = zug - zlg; yp = ylg - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_ts_rcv[count];

		zp = zug - zlg; yp = yug - ylg; xp = x - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_tn_rcv[count];

		count++;
	}
	//y edge
	count = 0;
	for(y = yl; y <= yu; y++)		
	{
		zp = zlg - zlg; yp = y - ylg; xp = xlg - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_bw_rcv[count];

		zp = zlg - zlg; yp = y - ylg; xp = xug - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_be_rcv[count];

		zp = zug - zlg; yp = y - ylg; xp = xlg - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_tw_rcv[count];

		zp = zug - zlg; yp = y - ylg; xp = xug - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_te_rcv[count];

		count++;
	}
	//z edge
	count = 0;
	for(z = zl; z <= zu; z++)
	{
		zp = z - zlg; yp = ylg - ylg; xp = xlg - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_sw_rcv[count];

		zp = z - zlg; yp = ylg - ylg; xp = xug - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_se_rcv[count];

		zp = z - zlg; yp = yug - ylg; xp = xlg - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_nw_rcv[count];

		zp = z - zlg; yp = yug - ylg; xp = xug - xlg;
		index = zp*NXYg + yp*NXg + xp;
		flag[index] = flag_ne_rcv[count];

		count++;
	}	
}
