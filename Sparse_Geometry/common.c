#include "common.h"
#include<stdio.h>
#include<stdlib.h>


int N, size;
//parameters used in physical field
double Lx, Ly, Lz;
double nu, U0, Db;
double Fx, Fy, Fz;
double Pin, Pout;
//dimensionless parameters
double Re = 0.1;
//parameters used in LBE simulation
double tau_f, wf, ci, rcc;
int e[Q][Dim] = 
{
	{0, 0, 0}, //0

	{1,  0, 0}, //1
	{-1, 0, 0}, //2
	{0,  1, 0}, //3
	{0, -1, 0}, //4
	{0,  0, 1}, //5
	{0,  0,-1}, //6

	{0,  1,  1},//7
	{0, -1, -1},//8
	{0, -1,  1},//9
	{0,  1, -1},//10
	{-1, 0, -1},//11
	{ 1, 0,  1},//12
	{-1, 0,  1},//13
	{ 1, 0, -1},//14
	{-1, 1,  0},//15
	{ 1,-1,  0},//16
	{-1, -1, 0},//17
	{ 1,  1, 0}//18
};
//parameters used in computation
unsigned int  n = 0, nmax;
double dx, dt;
double gama, rgama;
double sum_u_o = 0.0, sum_v_o = 0.0, sum_w_o = 0.0;

// host address
double *f, *fold, *fnew;
double *sf;
int  *node_index;
// device address
double *f_dev;  
double *F_dev;
int *node_dev;
double *sf_dev;
