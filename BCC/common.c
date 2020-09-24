#include "common.h"
#include<stdio.h>
#include<stdlib.h>

//parameters used in physical field
double nu, U0, Db;
double Fx, Fy, Fz;
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
int  n = 0, nmax;
double dx, dt;
double gama, rgama;
double sum_u_o = 0.0;

// CPUÖÐÄÚ´æ·ÖÅä
double *f, *fold, *fnew;
double *sf;
double *nc;
int   *flag;
// device adress
double *f_dev;  
double *F_dev;
int   *flag_dev;
double *sf_dev;
double *nc_dev;
