#ifndef __COMMON_H_
#define __COMMON_H_

#define  NX   90   //Grid in x- direciton
#define  NY   90   //Grid in y- direciton
#define  NZ   90   //Grid in z-direction
#define  NX1 (NX+1)  //Number of the grid in x-direction
#define  NY1 (NY+1)  //Number of the grid in y-direction
#define  NZ1 (NZ+1)  //Number of the grid in z-direction
#define  Lx 1.0      //Length of the domain in x-direction
#define  Ly (Lx*NY/NX)  //Length of the domain in y-direction
#define  Lz (Lx*NZ/NX)  //Length of the domain in z-direction
//------------------------------------------------------------------------------------------GPU
#define N16     16
#define NX2     ((NX1/16+1)*16+N16)
#define NY2     (NY+3)
#define NZ2     (NZ+3)
#define size    (NZ2*NY2*NX2)
#define NR      9
#define Dim     3
#define GID(z,y,x) (z*NX2*NY2+y*NX2+x)
////////////////////////////////////////////////////////////////////////////////////////////////////
#define BX      128 
#define BY      1    
#define BCX     64
////////////////////////////////////////////////////////////////////////////////////////////////////
#define Q   19
#define TP  100
#define PI  (4.0*atan(1.0)) 

//parameters used in physical field
extern double nu, U0, Db;
extern double Fx, Fy, Fz;
//dimensionless parameters
extern double Re;
//parameters used in LBE simulation
extern double tau_f, wf, ci, rcc;
extern double gama, rgama;
extern int e[Q][Dim];
//parameters used in computation
extern int  n, nmax;
extern double dx, dt;
extern double sum_u_o;

/////////////////////////////////////////////////////////////////////////////////////////////////
extern double *f, *fold, *fnew;
extern double *sf;
extern double *nc;
extern int   *flag;
//---------------------------------------------------------------------------------------------
// device adress
extern double *f_dev;  
extern double *F_dev;
extern int   *flag_dev;
extern double *sf_dev;
extern double *nc_dev;
//-----------------------------------------------------------------------------------------------------------------
#endif
