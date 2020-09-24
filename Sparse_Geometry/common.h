#ifndef __COMMON_H_
#define __COMMON_H_


//------------------------------------------------------------------------------------------GPU
#define Dim     3
////////////////////////////////////////////////////////////////////////////////////////////////////
#define BX      128
////////////////////////////////////////////////////////////////////////////////////////////////////
#define Q   19
#define TP  1000
#define PI  (4.0*atan(1.0)) 

//parameters used in physical field
extern int N, size;
extern double Lx, Ly, Lz;
extern double nu, U0, Db;
extern double Fx, Fy, Fz;
extern double Pin, Pout;
//dimensionless parameters
extern double Re;
//parameters used in LBE simulation
extern double tau_f, wf, ci, rcc;
extern double gama, rgama;
extern int e[Q][Dim];
//parameters used in computation
extern unsigned int  n, nmax;
extern double dx, dt;
extern double sum_u_o, sum_v_o, sum_w_o;

/////////////////////////////////////////////////////////////////////////////////////////////////
extern double *f, *fold, *fnew;
extern double *sf;
extern int  *node_index;
//---------------------------------------------------------------------------------------------
// device address
extern double *f_dev;  
extern double *F_dev;
extern int *node_dev;
extern double *sf_dev;
//-----------------------------------------------------------------------------------------------------------------
#endif
