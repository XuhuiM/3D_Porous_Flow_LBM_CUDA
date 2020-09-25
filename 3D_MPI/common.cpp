#include "common.h"
#include <cmath>

int iStep, MaxStep, tStat;

//domain size
int NX, NY, NZ, NXg, NYg, NZg;
int NXY, NXYZ;
int NXYg, NXYZg;
int xmin, xmax, xl, xu, xlg, xug;
int ymin, ymax, yl, yu, ylg, yug;
int zmin, zmax, zl, zu, zlg, zug;
int GNX, GNY, GNZ, GNXYZ;
int GNXY;
int xsize, ysize, zsize;
int xedge, yedge, zedge;
int xDVsize, yDVsize, zDVsize;
int xedgeDV, yedgeDV, zedgeDV;

//fluid parameters
double R;
double Lx, Ly, Lz;
double Re, U0, nu;
double L, dx, dt, tau, omega;
double *rho, *u, *v, *w;
double *rho_o, *u_o, *v_o, *w_o;
int    *flag;
double uConv;
double Fx;

//LBMParameters
int e[Q][ndim] = {
    { 0,  0,  0},

    { 1,  0,  0},
    {-1,  0,  0},
    { 0,  1,  0},
    { 0, -1,  0},
    { 0,  0,  1},
    { 0,  0, -1},

    { 0,  1,  1},
    { 0, -1, -1},
    { 0, -1,  1},
    { 0,  1, -1},

    {-1,  0, -1},
    { 1,  0,  1},
    {-1,  0,  1},
    { 1,  0, -1},

    {-1,  1,  0},
    { 1, -1,  0},
    {-1, -1,  0},
    { 1,  1,  0}
};
double Wi[Q] = {
	1.0/3, 
    	1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18, 1.0/18,
    	1.0/36, 1.0/36, 1.0/36, 1.0/36,
    	1.0/36, 1.0/36, 1.0/36, 1.0/36, 
    	1.0/36, 1.0/36, 1.0/36, 1.0/36
};
int re[Q] = {
	0,
    	2, 1, 4, 3, 6, 5, 
    	8, 7, 10, 9,
    	12, 11, 14, 13,
    	16, 15, 18, 17
};
double rcc = 3.0;
double *f, *fp;

//Communication parameters
int nproc, proc, vproc;
MPI_Comm MPI_COMM_VGRID;
int mpi_xdim, mpi_ydim, mpi_zdim;
int east, west, north, south, top, bottom;
int ne, nw, se, sw, te, tw, be, bw, tn, ts, bn, bs;
int master = 0;
int mpi_dim = 3;

// x face
double *f_west_snd, *f_east_snd, *f_west_rcv, *f_east_rcv;
int    *flag_west_snd, *flag_east_snd, *flag_west_rcv, *flag_east_rcv; 
// y face
double *f_south_snd, *f_north_snd, *f_south_rcv, *f_north_rcv;
int    *flag_south_snd, *flag_north_snd, *flag_south_rcv, *flag_north_rcv;
// z face
double *f_bot_snd, *f_top_snd, *f_bot_rcv, *f_top_rcv;
int    *flag_bot_snd, *flag_top_snd, *flag_bot_rcv, *flag_top_rcv;

// x edge
double *f_tn_snd, *f_tn_rcv;
double *f_bn_snd, *f_bn_rcv;
double *f_ts_snd, *f_ts_rcv;
double *f_bs_snd, *f_bs_rcv;
int    *flag_tn_snd, *flag_tn_rcv;
int    *flag_bn_snd, *flag_bn_rcv;
int    *flag_ts_snd, *flag_ts_rcv;
int    *flag_bs_snd, *flag_bs_rcv;

// y edge
double *f_te_snd, *f_te_rcv;
double *f_be_snd, *f_be_rcv;
double *f_tw_snd, *f_tw_rcv;
double *f_bw_snd, *f_bw_rcv;
int    *flag_te_snd, *flag_te_rcv;
int    *flag_be_snd, *flag_be_rcv;
int    *flag_tw_snd, *flag_tw_rcv;
int    *flag_bw_snd, *flag_bw_rcv; 

// z edge
double *f_ne_snd, *f_ne_rcv;
double *f_se_snd, *f_se_rcv;
double *f_nw_snd, *f_nw_rcv;
double *f_sw_snd, *f_sw_rcv;
int    *flag_ne_snd, *flag_ne_rcv;
int    *flag_se_snd, *flag_se_rcv;
int    *flag_nw_snd, *flag_nw_rcv;
int    *flag_sw_snd, *flag_sw_rcv;

//equilibrium functions
double feq(double RHO, double U, double V, double W, int k)
{
	double FEQ, euvw, uvw;
	euvw = e[k][0]*U + e[k][1]*V + e[k][2]*W;
	uvw = U*U + V*V + W*W;
	FEQ = Wi[k]*(RHO + 3.0*euvw + 4.5*euvw*euvw - 1.5*uvw);
	return FEQ;
}
