/////// allocate or delete memory for the pointers////////////////
#include "common.h"

void memalloc(int FLAG)
{
	
	//malloc arrarys
	if(FLAG == 1)
	{
		NX = xu - xl + 1;
		NY = yu - yl + 1;
		NZ = zu - zl + 1;
		NXg = NX + 2;
		NYg = NY + 2;
		NZg = NZ + 2;
		NXYg = NXg*NYg;
		NXYZg = NXg*NYg*NZg;
	
		//Face
		xsize  = NZ*NY;
    		xDVsize = xsize*Q;
    		ysize  = NZ*NX;
    		yDVsize = ysize*Q;
    		zsize  = NY*NX;
    		zDVsize = zsize*Q;

    		//Edge
		xedge = NX;
		yedge = NY;
		zedge = NZ;
    		xedgeDV = NX*Q;
    		yedgeDV = NY*Q;
    		zedgeDV = NZ*Q;

		GNX = xmax - xmin + 1;
		GNY = ymax - ymin + 1;
		GNZ = zmax - zmin + 1;
		GNXY = GNX*GNY;
		GNXYZ = GNX*GNY*GNZ;
		
		rho = new double[NXYZg];
		u   = new double[NXYZg];
		v   = new double[NXYZg];
		w   = new double[NXYZg];
		flag = new int[NXYZg];

		f =  new double[Q*NXYZg];
		fp = new double[Q*NXYZg];

		rho_o = new double[GNXYZ];
		u_o   = new double[GNXYZ];
		v_o   = new double[GNXYZ];
		w_o   = new double[GNXYZ];

		// Allocation for exchange of data in x direction
		f_west_snd = new double[xDVsize];
       		f_east_snd = new double[xDVsize];
       		f_west_rcv = new double[xDVsize];
       		f_east_rcv = new double[xDVsize];
		flag_west_snd = new int[xsize];
       		flag_east_snd = new int[xsize];
       		flag_west_rcv = new int[xsize];
       		flag_east_rcv = new int[xsize];

		// Allocation for exchange of data in y direction
       		f_south_snd = new double[yDVsize];
       		f_north_snd = new double[yDVsize];
       		f_south_rcv = new double[yDVsize];
       		f_north_rcv = new double[yDVsize];
       		flag_south_snd = new int[ysize];
       		flag_north_snd = new int[ysize];
       		flag_south_rcv = new int[ysize];
       		flag_north_rcv = new int[ysize];

       		// Allocation for exchange of data in the z direction
       		f_top_snd = new double[zDVsize];
       		f_bot_snd = new double[zDVsize];
       		f_top_rcv = new double[zDVsize];
       		f_bot_rcv = new double[zDVsize];
       		flag_top_snd = new int[zsize];
       		flag_bot_snd = new int[zsize];
       		flag_top_rcv = new int[zsize];
       		flag_bot_rcv = new int[zsize];

		// Allocation for data exchange in edges along the x direction
       		f_tn_snd = new double[xedgeDV];
       		f_ts_snd = new double[xedgeDV];
       		f_bn_snd = new double[xedgeDV];
       		f_bs_snd = new double[xedgeDV];
       		f_tn_rcv = new double[xedgeDV];
       		f_ts_rcv = new double[xedgeDV];
       		f_bn_rcv = new double[xedgeDV];
       		f_bs_rcv = new double[xedgeDV];
       		flag_tn_snd = new int[xedge];
       		flag_ts_snd = new int[xedge];
       		flag_bn_snd = new int[xedge];
       		flag_bs_snd = new int[xedge];
       		flag_tn_rcv = new int[xedge];
       		flag_ts_rcv = new int[xedge];
       		flag_bn_rcv = new int[xedge];
       		flag_bs_rcv = new int[xedge];

      		// Allocation for data exchange in edges along the y direction
       		f_te_snd = new double[yedgeDV];
       		f_tw_snd = new double[yedgeDV];
       		f_be_snd = new double[yedgeDV];
       		f_bw_snd = new double[yedgeDV];
       		f_te_rcv = new double[yedgeDV];
       		f_tw_rcv = new double[yedgeDV];
       		f_be_rcv = new double[yedgeDV];
       		f_bw_rcv = new double[yedgeDV];
       		flag_te_snd = new int[yedge];
       		flag_tw_snd = new int[yedge];
       		flag_be_snd = new int[yedge];
       		flag_bw_snd = new int[yedge];
       		flag_te_rcv = new int[yedge];
       		flag_tw_rcv = new int[yedge];
       		flag_be_rcv = new int[yedge];
       		flag_bw_rcv = new int[yedge];

       		// Allocation for data exchange in edges along the z direction
       		f_ne_snd = new double[zedgeDV];
       		f_nw_snd = new double[zedgeDV];
       		f_se_snd = new double[zedgeDV];
       		f_sw_snd = new double[zedgeDV];
       		f_ne_rcv = new double[zedgeDV];
       		f_nw_rcv = new double[zedgeDV];
       		f_se_rcv = new double[zedgeDV];
       		f_sw_rcv = new double[zedgeDV];
       		flag_ne_snd = new int[zedge];
       		flag_nw_snd = new int[zedge];
       		flag_se_snd = new int[zedge];
       		flag_sw_snd = new int[zedge];
       		flag_ne_rcv = new int[zedge];
       		flag_nw_rcv = new int[zedge];
       		flag_se_rcv = new int[zedge];
       		flag_sw_rcv = new int[zedge];

	}
	else
	{
		delete [] rho;
        	delete [] u;
        	delete [] v;
        	delete [] w;
		delete [] flag;

		delete [] f;
		delete [] fp;

		delete [] rho_o;
        	delete [] u_o;
        	delete [] v_o;
        	delete [] w_o;
				
		// free data in the x direction
        	delete [] f_west_snd;
        	delete [] f_east_snd;
        	delete [] f_west_rcv;
        	delete [] f_east_rcv;
        	delete [] flag_west_snd;
        	delete [] flag_east_snd;
        	delete [] flag_west_rcv;
        	delete [] flag_east_rcv;

        	// free data in the y direction
        	delete [] f_south_snd;
        	delete [] f_north_snd;
        	delete [] f_south_rcv;
        	delete [] f_north_rcv;
        	delete [] flag_south_snd;
        	delete [] flag_north_snd;
        	delete [] flag_south_rcv;
        	delete [] flag_north_rcv;

        	// free data in the z direction
        	delete [] f_top_snd;
        	delete [] f_bot_snd;
        	delete [] f_top_rcv;
        	delete [] f_bot_rcv;
        	delete [] flag_top_snd;
        	delete [] flag_bot_snd;
        	delete [] flag_top_rcv;
        	delete [] flag_bot_rcv;

        	// free data in edges along the x direction
        	delete [] f_tn_snd;
        	delete [] f_ts_snd;
        	delete [] f_bn_snd;
        	delete [] f_bs_snd;
        	delete [] f_tn_rcv;
        	delete [] f_ts_rcv;
        	delete [] f_bn_rcv;
        	delete [] f_bs_rcv;
        	delete [] flag_tn_snd;
        	delete [] flag_ts_snd;
        	delete [] flag_bn_snd;
        	delete [] flag_bs_snd;
        	delete [] flag_tn_rcv;
        	delete [] flag_ts_rcv;
        	delete [] flag_bn_rcv;
        	delete [] flag_bs_rcv;

        	// free data in edges along the y direction
        	delete [] f_te_snd;
        	delete [] f_tw_snd;
        	delete [] f_be_snd;
        	delete [] f_bw_snd;
        	delete [] f_te_rcv;
        	delete [] f_tw_rcv;
        	delete [] f_be_rcv;
        	delete [] f_bw_rcv;
        	delete [] flag_te_snd;
        	delete [] flag_tw_snd;
        	delete [] flag_be_snd;
        	delete [] flag_bw_snd;
        	delete [] flag_te_rcv;
        	delete [] flag_tw_rcv;
        	delete [] flag_be_rcv;
        	delete [] flag_bw_rcv;

        	// free data in edges along the z direction
        	delete [] f_ne_snd;
        	delete [] f_nw_snd;
        	delete [] f_se_snd;
        	delete [] f_sw_snd;
        	delete [] f_ne_rcv;
        	delete [] f_nw_rcv;
        	delete [] f_se_rcv;
        	delete [] f_sw_rcv;
        	delete [] flag_ne_snd;
        	delete [] flag_nw_snd;
        	delete [] flag_se_snd;
        	delete [] flag_sw_snd;
        	delete [] flag_ne_rcv;
        	delete [] flag_nw_rcv;
        	delete [] flag_se_rcv;
        	delete [] flag_sw_rcv;
	}
}
