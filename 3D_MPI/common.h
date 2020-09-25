/// global variables used //////////////////////////////// 
#ifndef __COMMON_H_
#define __COMMON_H_

#include <cmath>
#include <sstream>
#include <string>
#include "mpi.h"

#define ndim 3
#define Q 19

extern int iStep, MaxStep, tStat;

//domain size
extern int NX, NY, NZ, NXg, NYg, NZg;
extern int NXY, NXYZ;
extern int NXYg, NXYZg;
extern int xmin, xmax, xl, xu, xlg, xug;
extern int ymin, ymax, yl, yu, ylg, yug;
extern int zmin, zmax, zl, zu, zlg, zug;
extern int GNX, GNY, GNZ, GNXYZ;
extern int GNXY;
extern int xsize, ysize, zsize;
extern int xedge, yedge, zedge;
extern int xDVsize, yDVsize, zDVsize;
extern int xedgeDV, yedgeDV, zedgeDV;

//fluid parameters
extern double R;
extern double Lx, Ly, Lz;
extern double Re, U0, nu;
extern double L, dx, dt, tau, omega;
extern double *rho, *u, *v, *w;
extern double *rho_o, *u_o, *v_o, *w_o;
extern int    *flag;
extern double Fx;
extern double uConv;

//LBMParameters
extern int e[Q][ndim];
extern double Wi[Q];
extern int re[Q];
extern double rcc;
extern double *f, *fp;

//Communication parameters
extern int nproc, proc, vproc;
extern MPI_Comm MPI_COMM_VGRID;
extern int mpi_xdim, mpi_ydim, mpi_zdim;
extern int east, west, north, south, top, bottom;
extern int ne, nw, se, sw, te, tw, be, bw, tn, ts, bn, bs;
extern int master;
extern int mpi_dim;

// x face
extern double *f_west_snd, *f_east_snd, *f_west_rcv, *f_east_rcv;
extern int    *flag_west_snd, *flag_east_snd, *flag_west_rcv, *flag_east_rcv; 
// y face
extern double *f_south_snd, *f_north_snd, *f_south_rcv, *f_north_rcv;
extern int    *flag_south_snd, *flag_north_snd, *flag_south_rcv, *flag_north_rcv;
// z face
extern double *f_bot_snd, *f_top_snd, *f_bot_rcv, *f_top_rcv;
extern int    *flag_bot_snd, *flag_top_snd, *flag_bot_rcv, *flag_top_rcv;

// x edge
extern double *f_tn_snd, *f_tn_rcv;
extern double *f_bn_snd, *f_bn_rcv;
extern double *f_ts_snd, *f_ts_rcv;
extern double *f_bs_snd, *f_bs_rcv;
extern int    *flag_tn_snd, *flag_tn_rcv;
extern int    *flag_bn_snd, *flag_bn_rcv;
extern int    *flag_ts_snd, *flag_ts_rcv;
extern int    *flag_bs_snd, *flag_bs_rcv;

// y edge
extern double *f_te_snd, *f_te_rcv;
extern double *f_be_snd, *f_be_rcv;
extern double *f_tw_snd, *f_tw_rcv;
extern double *f_bw_snd, *f_bw_rcv;
extern int    *flag_te_snd, *flag_te_rcv;
extern int    *flag_be_snd, *flag_be_rcv;
extern int    *flag_tw_snd, *flag_tw_rcv;
extern int    *flag_bw_snd, *flag_bw_rcv; 

// z edge
extern double *f_ne_snd, *f_ne_rcv;
extern double *f_se_snd, *f_se_rcv;
extern double *f_nw_snd, *f_nw_rcv;
extern double *f_sw_snd, *f_sw_rcv;
extern int    *flag_ne_snd, *flag_ne_rcv;
extern int    *flag_se_snd, *flag_se_rcv;
extern int    *flag_nw_snd, *flag_nw_rcv;
extern int    *flag_sw_snd, *flag_sw_rcv;

// void funcs
void parameters();
void geo();
void vgrid();
void memalloc(int FLAG);
void init();
void collision();
void postcollision();
void flagex();
void relaxstats();
void stream();
void boundary();
void datasave(int step);
double feq(double RHO, double U, double V, double W, int k);

#endif
