#include "lb.h"
#include "common.h"

__constant__ int e_d[Q][Dim] = 
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

__constant__ int re_d[Q] = 
{
	0,
	2,
	1,
	4,
	3,
	6,
	5,


	8,
	7,
	10,
	9,
	12,
	11,
	14,
	13,
	16,
	15,
	18,
	17
};

__constant__ double w_d[Q] = 
{
	1.0/3,

	1.0/18,
	1.0/18,
	1.0/18,
	1.0/18,
	1.0/18,
	1.0/18,

	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36,
	1.0/36
};

__constant__ double r0 = 1.0/19;
__constant__ double r1 = 5.0/399;
__constant__ double r2 = 1.0/21;
__constant__ double r3 = 11.0/2394;
__constant__ double r4 = 1.0/63;
__constant__ double r5 = 1.0/10;
__constant__ double r6 = 1.0/18;
__constant__ double r7 = 1.0/36;
__constant__ double r8 = 4.0/1197;
__constant__ double r9 = 1.0/252;
__constant__ double r10 = 1.0/40;
__constant__ double r11 = 1.0/12;
__constant__ double r12 = 1.0/4;
__constant__ double r13 = 1.0/72;
__constant__ double r14 = 1.0/24;
__constant__ double r15 = 1.0/8;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

__global__ void Evol_flow(double rgama, double *sf_d, double dt, double Fx, double Fy, double Fz, double *f_d, double *F_d)
{
	double P, U, V, W, UVW;
	int tx;
	int bx, by, bz;
	int z, y, x, k;	
	
	double mf[Q];
	double f0[BX];
	double f1[BX], f2[BX], f3[BX], f4[BX], f5[BX], f6[BX];
	double f7[BX], f8[BX], f9[BX], f10[BX], f11[BX], f12[BX], f13[BX], f14[BX], f15[BX], f16[BX], f17[BX], f18[BX];
	__shared__ double F0[BX];
	__shared__ double F1[BX], F2[BX], F3[BX], F4[BX], F5[BX], F6[BX];
	__shared__ double F7[BX], F8[BX], F9[BX], F10[BX], F11[BX], F12[BX], F13[BX], F14[BX], F15[BX], F16[BX], F17[BX], F18[BX];

	tx = threadIdx.x;
    bx = blockIdx.x;
	by = blockIdx.y;
	bz = blockIdx.z;	
	x = N16+bx*BX+tx;
	y = 1+by;
    z = 1+bz;	
	k = GID(z,y,x);


	if(x <= N16+NX1)
	{

		f0[tx]  = f_d[k+0*size]; 
		f1[tx]  = f_d[k+1*size]; 
		f2[tx]  = f_d[k+2*size]; 
		f3[tx]  = f_d[k+3*size]; 
		f4[tx]  = f_d[k+4*size]; 
		f5[tx]  = f_d[k+5*size]; 
		f6[tx]  = f_d[k+6*size]; 
		f7[tx]  = f_d[k+7*size]; 
		f8[tx]  = f_d[k+8*size]; 
		f9[tx]  = f_d[k+9*size]; 
		f10[tx] = f_d[k+10*size]; 
		f11[tx] = f_d[k+11*size]; 
		f12[tx] = f_d[k+12*size]; 
		f13[tx] = f_d[k+13*size]; 
		f14[tx] = f_d[k+14*size]; 
		f15[tx] = f_d[k+15*size]; 
		f16[tx] = f_d[k+16*size]; 
		f17[tx] = f_d[k+17*size]; 
		f18[tx] = f_d[k+18*size]; 
		
		//f-mf///////////////////////////
		mf[0]  = f0[tx] + f1[tx] + f2[tx] + f3[tx] + f4[tx] + f5[tx] + f6[tx] + f7[tx] + f8[tx] + f9[tx] + f10[tx] + f11[tx] + f12[tx] + f13[tx] + f14[tx] + f15[tx] + f16[tx] + f17[tx] + f18[tx]; 
 		mf[1]  = -30.0*f0[tx] - 11.0*f1[tx] - 11.0*f2[tx] - 11.0*f3[tx] - 11.0*f4[tx] - 11.0*f5[tx]
		       	-11.0*f6[tx] + 8.0*f7[tx] + 8.0*f8[tx] + 8.0*f9[tx] + 8.0*f10[tx] + 8.0*f11[tx] + 8.0*f12[tx] + 8.0*f13[tx] + 8.0*f14[tx] + 8.0*f15[tx] + 8.0*f16[tx] + 8.0*f17[tx] + 8.0*f18[tx];
                mf[2]  = 12.0*f0[tx] - 4.0*f1[tx] - 4.0*f2[tx] - 4.0*f3[tx] - 4.0*f4[tx] - 4.0*f5[tx] - 4.0*f6[tx] + f7[tx] + f8[tx] + f9[tx] 
			+ f10[tx] + f11[tx] + f12[tx] + f13[tx] + f14[tx] + f15[tx] + f16[tx] + f17[tx] + f18[tx];
              	mf[3]  = f1[tx] - f2[tx] - f11[tx] + f12[tx] - f13[tx] + f14[tx] - f15[tx] + f16[tx] - f17[tx] + f18[tx];
                mf[4]  = -4.0*f1[tx] + 4.0*f2[tx]  - f11[tx] + f12[tx] - f13[tx] + f14[tx] - f15[tx] + f16[tx] - f17[tx] + f18[tx];
                mf[5]  = f3[tx] - f4[tx] + f7[tx] - f8[tx] - f9[tx] + f10[tx] + f15[tx] - f16[tx] - f17[tx] + f18[tx];
                mf[6]  = -4.0*f3[tx] + 4.0*f4[tx] + f7[tx] - f8[tx] - f9[tx] + f10[tx] + f15[tx] - f16[tx] - f17[tx] + f18[tx];
                mf[7]  = f5[tx] - f6[tx] + f7[tx] - f8[tx] + f9[tx] + f12[tx] - f11[tx] - f10[tx] + f13[tx] - f14[tx];
                mf[8]  = f12[tx] - f11[tx] - f10[tx] + f13[tx] - f14[tx] - 4.0*f5[tx] + 4.0*f6[tx] + f7[tx] - f8[tx] + f9[tx];
                mf[9]  = 2.0*f1[tx] - 2.0*f10[tx] + f11[tx] + f12[tx] + f13[tx] + f14[tx] + f15[tx] + f16[tx] + f17[tx] + f18[tx] + 2.0*f2[tx] - f3[tx] - f4[tx] - f5[tx] - f6[tx] - 2.0*f7[tx] - 2.0*f8[tx] - 2.0*f9[tx];
                mf[10] = f11[tx] - 2.0*f10[tx] - 4.0*f1[tx] + f12[tx] + f13[tx] + f14[tx]
		       	+ f15[tx] + f16[tx] + f17[tx] + f18[tx] - 4.0*f2[tx] + 2.0*f3[tx] + 2.0*f4[tx] + 2.0*f5[tx] + 2.0*f6[tx] - 2.0*f7[tx] - 2.0*f8[tx] - 2.0*f9[tx];
                mf[11] = f15[tx] - f12[tx] - f13[tx] - f14[tx] - f11[tx] + f16[tx] + f17[tx] + f18[tx] + f3[tx] + f4[tx] - f5[tx] - f6[tx];
                mf[12] = f15[tx] - f12[tx] - f13[tx] - f14[tx] - f11[tx] + f16[tx] + f17[tx] + f18[tx] - 2.0*f3[tx] - 2.0*f4[tx] + 2.0*f5[tx] + 2.0*f6[tx];
                mf[13] = f17[tx] - f16[tx] - f15[tx] + f18[tx];
                mf[14] = f7[tx] - f10[tx] + f8[tx] - f9[tx];
                mf[15] = f11[tx] + f12[tx] - f13[tx] - f14[tx];
                mf[16] = f11[tx] - f12[tx] + f13[tx] - f14[tx] - f15[tx] + f16[tx] - f17[tx] + f18[tx];
                mf[17] = f10[tx] - f15[tx] + f16[tx] + f17[tx] - f18[tx] + f7[tx] - f8[tx] - f9[tx];
                mf[18] = f10[tx] - f11[tx] + f12[tx] + f13[tx] - f14[tx] - f7[tx] + f8[tx] - f9[tx];

		//macroscopic variables/////////////////////////////////////////////////////////////////////////////
		U   = f1[tx] + f12[tx] + f14[tx] + f16[tx] + f18[tx] - f2[tx] - f11[tx] - f13[tx] - f15[tx] - f17[tx] + 0.5*dt*Fx;
		V   = f3[tx] + f7[tx]  + f10[tx] + f15[tx] + f18[tx] - f4[tx] - f8[tx]  - f9[tx]  - f16[tx] - f17[tx] + 0.5*dt*Fy;
		W   = f5[tx] + f7[tx]  + f9[tx]  + f12[tx] + f13[tx] - f6[tx] - f8[tx]  - f10[tx] - f11[tx] - f14[tx] + 0.5*dt*Fz;
		UVW = U*U + V*V + W*W;
		P = (f0[tx] + f1[tx] + f2[tx] + f3[tx] + f4[tx] + f5[tx] + f6[tx] + f7[tx] + f8[tx] + f9[tx] + f10[tx] + f11[tx] + f12[tx] + f13[tx] 
				+ f14[tx] + f15[tx] + f16[tx] + f17[tx] + f18[tx]);
		
		//collision//-------------------------------------------------------------------------------------------------------------------------------
		mf[0]  = (mf[0] - sf_d[0]*( mf[0] - MEQ_0(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[0])*F_0(U,V,W,Fx,Fy,Fz,rgama);
		mf[1]  = (mf[1] - sf_d[1]*( mf[1] - MEQ_1(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[1])*F_1(U,V,W,Fx,Fy,Fz,rgama);
		mf[2]  = (mf[2] - sf_d[2]*( mf[2] - MEQ_2(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[2])*F_2(U,V,W,Fx,Fy,Fz,rgama);
		mf[3]  = (mf[3] - sf_d[3]*( mf[3] - MEQ_3(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[3])*F_3(U,V,W,Fx,Fy,Fz,rgama);
		mf[4]  = (mf[4] - sf_d[4]*( mf[4] - MEQ_4(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[4])*F_4(U,V,W,Fx,Fy,Fz,rgama);
		mf[5]  = (mf[5] - sf_d[5]*( mf[5] - MEQ_5(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[5])*F_5(U,V,W,Fx,Fy,Fz,rgama);
		mf[6]  = (mf[6] - sf_d[6]*( mf[6] - MEQ_6(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[6])*F_6(U,V,W,Fx,Fy,Fz,rgama);
		mf[7]  = (mf[7] - sf_d[7]*( mf[7] - MEQ_7(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[7])*F_7(U,V,W,Fx,Fy,Fz,rgama);
		mf[8]  = (mf[8] - sf_d[8]*( mf[8] - MEQ_8(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[8])*F_8(U,V,W,Fx,Fy,Fz,rgama);
		mf[9]  = (mf[9] - sf_d[9]*( mf[9] - MEQ_9(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[9])*F_9(U,V,W,Fx,Fy,Fz,rgama);
		mf[10] = (mf[10]-sf_d[10]*(mf[10] -MEQ_10(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[10])*F_10(U,V,W,Fx,Fy,Fz,rgama);
		mf[11] = (mf[11]-sf_d[11]*(mf[11] -MEQ_11(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[11])*F_11(U,V,W,Fx,Fy,Fz,rgama);
		mf[12] = (mf[12]-sf_d[12]*(mf[12] -MEQ_12(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[12])*F_12(U,V,W,Fx,Fy,Fz,rgama);
		mf[13] = (mf[13]-sf_d[13]*(mf[13] -MEQ_13(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[13])*F_13(U,V,W,Fx,Fy,Fz,rgama);
		mf[14] = (mf[14]-sf_d[14]*(mf[14] -MEQ_14(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[14])*F_14(U,V,W,Fx,Fy,Fz,rgama);
		mf[15] = (mf[15]-sf_d[15]*(mf[15] -MEQ_15(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[15])*F_15(U,V,W,Fx,Fy,Fz,rgama);
		mf[16] = (mf[16]-sf_d[16]*(mf[16] -MEQ_16(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[16])*F_16(U,V,W,Fx,Fy,Fz,rgama);
		mf[17] = (mf[17]-sf_d[17]*(mf[17] -MEQ_17(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[17])*F_17(U,V,W,Fx,Fy,Fz,rgama);
		mf[18] = (mf[18]-sf_d[18]*(mf[18] -MEQ_18(P, U, V, W, rgama)))+dt*(1.0 -0.5*sf_d[18])*F_18(U,V,W,Fx,Fy,Fz,rgama);

		//--mf - f --//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		F0[tx]  = r0*mf[0] - r1*mf[1] + r2*mf[2];
                F1[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[3] - r5*mf[4] + r6*mf[9] - r6*mf[10];
                F2[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[3] + r5*mf[4] + r6*mf[9] - r6*mf[10]; 
                F3[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[5] - r5*mf[6] - r7*mf[9] + r7*mf[10] + r11*mf[11] - r11*mf[12]; 
                F4[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[5] + r5*mf[6] - r7*mf[9] + r7*mf[10] + r11*mf[11] - r11*mf[12];
                F5[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[7] - r5*mf[8] - r7*mf[9] + r7*mf[10] - r11*mf[11] + r11*mf[12];
                F6[tx]  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[7] + r5*mf[8] - r7*mf[9] + r7*mf[10] - r11*mf[11] + r11*mf[12];
                F7[tx]  = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[5] + r10*mf[6] + r5*mf[7] + r10*mf[8] - r6*mf[9]  - r7*mf[10] + r12*mf[14] + r15*mf[17] - r15*mf[18];
                F8[tx]  = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[5] - r10*mf[6] - r5*mf[7] - r10*mf[8] - r6*mf[9]  - r7*mf[10] + r12*mf[14] - r15*mf[17] + r15*mf[18];
                F9[tx]  = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[5] - r10*mf[6] + r5*mf[7] + r10*mf[8] - r6*mf[9]  - r7*mf[10] - r12*mf[14] - r15*mf[17] - r15*mf[18];
                F10[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[5] + r10*mf[6] - r5*mf[7] - r10*mf[8] - r6*mf[9]  - r7*mf[10] - r12*mf[14] + r15*mf[17] + r15*mf[18];
                F11[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] - r5*mf[7] - r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] + r12*mf[15] + r15*mf[16] - r15*mf[18];
                F12[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] + r5*mf[7] + r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] + r12*mf[15] - r15*mf[16] + r15*mf[18];
                F13[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] + r5*mf[7] + r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] - r12*mf[15] + r15*mf[16] + r15*mf[18]; 
                F14[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] - r5*mf[7] - r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] - r12*mf[15] - r15*mf[16] - r15*mf[18]; 
                F15[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] + r5*mf[5] + r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] - r12*mf[13] - r15*mf[16] - r15*mf[17]; 
                F16[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] - r5*mf[5] - r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] - r12*mf[13] + r15*mf[16] + r15*mf[17]; 
                F17[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] - r5*mf[5] - r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] + r12*mf[13] - r15*mf[16] + r15*mf[17];
                F18[tx] = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] + r5*mf[5] + r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] + r12*mf[13] + r15*mf[16] - r15*mf[17];
		
		__syncthreads();

		//  streaming  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		F_d[k           + 0*size] = F0[tx];


		F_d[k + NX2     + 3*size] = F3[tx];
		F_d[k - NX2     + 4*size] = F4[tx];
		F_d[k + NX2*NY2 + 5*size] = F5[tx];
		F_d[k - NX2*NY2 + 6*size] = F6[tx];


		F_d[k + NX2 + NX2*NY2 + 7*size] = F7[tx];
		F_d[k - NX2 - NX2*NY2 + 8*size] = F8[tx];
		F_d[k - NX2 + NX2*NY2 + 9*size] = F9[tx];
		F_d[k + NX2 - NX2*NY2 +10*size] = F10[tx];

		if(tx != 0)
		{
			F_d[k           +  1*size] =  F1[tx-1];
			F_d[k + NX2*NY2 + 12*size] = F12[tx-1];
			F_d[k - NX2*NY2 + 14*size] = F14[tx-1];
			F_d[k - NX2     + 16*size] = F16[tx-1];
			F_d[k + NX2     + 18*size] = F18[tx-1];
		}

		if(tx == BX-1)
		{
			F_d[k + 1           +  1*size] = F1[tx];
			F_d[k + 1 + NX2*NY2 + 12*size] = F12[tx];
			F_d[k + 1 - NX2*NY2 + 14*size] = F14[tx];
			F_d[k + 1 - NX2     + 16*size] = F16[tx];
			F_d[k + 1 + NX2     + 18*size] = F18[tx];
		}

		if(tx != BX-1)
		{	
			F_d[k           +  2*size] = F2[tx+1];
			F_d[k - NX2*NY2 + 11*size] = F11[tx+1];
			F_d[k + NX2*NY2 + 13*size] = F13[tx+1];
			F_d[k + NX2     + 15*size] = F15[tx+1];
			F_d[k - NX2     + 17*size] = F17[tx+1];
		}

		if(tx == 0)
		{
			F_d[k - 1           +  2*size] = F2[tx];
			F_d[k - 1 - NX2*NY2 + 11*size] = F11[tx];
			F_d[k - 1 + NX2*NY2 + 13*size] = F13[tx];
			F_d[k - 1 + NX2     + 15*size] = F15[tx];
			F_d[k - 1 - NX2     + 17*size] = F17[tx];        
		}
	}
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Bc_flow_BB(int *flag_d, double *F_d)
{
	int tx, bx, by, bz;
	int z, y, x, k;
	int zp, yp, xp, kp;
	int q;


	tx = threadIdx.x;
     	bx = blockIdx.x;
	by = blockIdx.y;
	bz = blockIdx.z;	
	x = N16+bx*BX+tx;
	y = 1+by;
       	z = 1+bz;	
	k = GID(z,y,x);

	if(x < N16+NX1)
	{
		if(flag_d[k] == 0)
		{
			for(q = 1; q < Q; q++)
			{
				xp = x - e_d[q][0]; yp = y - e_d[q][1]; zp = z - e_d[q][2];
				kp = GID(zp,yp,xp);
				if(flag_d[kp] == 1)
				{
					F_d[k + q*size] = F_d[kp + re_d[q]*size];
				}

			}
		}
	}

}

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Bc_flow_io(double *F_d)       
{
	int ty, bx, bz;
	int z, y, x, k;
	int zp, yp, xp, kp;

	ty = threadIdx.x; 
	bx = blockIdx.x;
	bz = blockIdx.y;
	y = 1 + bx*BCX + ty;
	z = 1 + bz;

	if(y <= NY1)
	{
		//inlet
		x = N16;
		k = GID(z, y, x);
		// c_1
		xp = N16+NX1;
		yp = y; zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 1*size] = F_d[kp + 1*size];

		// c_12
		yp = y;
       		if(z > 1) zp = z;
		else zp = NZ1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 12*size] = F_d[kp + 12*size];

		// c_14
		yp = y;
		if(z < NZ1) zp = z;
		else zp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 14*size] = F_d[kp + 14*size];

		// c_16
		if(y < NY1) yp = y;
		else yp = 0; 
		zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 16*size] = F_d[kp + 16*size];

		// c_18
		if(y > 1) yp = y;
		else yp = NY1+1; 
		zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 18*size] = F_d[kp + 18*size];

		//outlet
		x = N16+NX;
		k = GID(z, y, x);
		// c_2
		xp = N16-1;
		yp = y; zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 2*size] = F_d[kp + 2*size];

		// c_11
		yp = y;
       		if(z < NZ1) zp = z;
		else zp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 11*size] = F_d[kp + 11*size];

		// c_13
		yp = y;
		if(z > 1) zp = z;
		else zp = NZ1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 13*size] = F_d[kp + 13*size];

		// c_15
		if(y > 1) yp = y;
		else yp = NY1+1; 
		zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 15*size] = F_d[kp + 15*size];

		// c_17
		if(y < NY1) yp = y;
		else yp = 0; 
		zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 17*size] = F_d[kp + 17*size];	
	}
}
//------------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Bc_flow_fb(double *F_d)       
{
	int tx, bx, bz;
	int z, y, x, k;
	int zp, yp, xp, kp;

	tx = threadIdx.x; 
	bx = blockIdx.x;
	bz = blockIdx.y;
	x = N16 + bx*BCX + tx;
	z = 1 + bz;

	if(x < N16+NX1)
	{
		//front
		y = 1;
		k = GID(z, y, x);
		// c_3
		yp = NY1+1;
		xp = x; zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 3*size] = F_d[kp + 3*size];

		// c_7
		xp = x;
      		if(z > 1) zp = z;
		else zp = NZ1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 7*size] = F_d[kp + 7*size];

		// c_10
		xp = x;
		if(z < NZ1) zp = z;
		else zp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 10*size] = F_d[kp + 10*size];

		// c_15
		zp = z;
		if(x < N16+NX) xp = x;
		else xp = N16-1;
		kp = GID(zp, yp, xp);
		F_d[k + 15*size] = F_d[kp + 15*size];

		// c_18
		zp = z;
		if(x > N16) xp = x;
		else xp = N16+NX1;
		kp = GID(zp, yp, xp);
		F_d[k + 18*size] = F_d[kp + 18*size];

		//back
		y = NY1;
		k = GID(z, y, x);
		// c_4
		yp = 0;
		xp = x; zp = z;
		kp = GID(zp, yp, xp);
		F_d[k + 4*size] = F_d[kp + 4*size];

		// c_8
		xp = x;
    		if(z < NZ1) zp = z;
		else zp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 8*size] = F_d[kp + 8*size];

		// c_9
		xp = x;
		if(z > 1) zp = z;
		else zp = NZ1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 9*size] = F_d[kp + 9*size];

		// c_16
		zp = z;
		if(x > N16) xp = x;
		else xp = N16+NX1;
		kp = GID(zp, yp, xp);
		F_d[k + 16*size] = F_d[kp + 16*size];

		// c_17
		zp = z;
		if(x < N16+NX) xp = x;
		else xp = N16-1;
		kp = GID(zp, yp, xp);
		F_d[k + 17*size] = F_d[kp + 17*size];
	}
	

}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
__global__ void Bc_flow_ub(double *F_d)       
{
	int tx, bx, by;
	int z, y, x, k;
	int zp, yp, xp, kp;

	tx = threadIdx.x; 
	bx = blockIdx.x;
	by = blockIdx.y;
	x = N16 + bx*BCX + tx;
	y = 1 + by;

	if(x < N16+NX1)
	{
		//bottom
		z = 1;
		k = GID(z, y, x);
		// c_5
		zp = NZ1+1;
		xp = x; yp = y;
		kp = GID(zp, yp, xp);
		F_d[k + 5*size] = F_d[kp + 5*size];

		// c_7
		xp = x; 
		if(y > 1) yp = y;
		else yp = NY1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 7*size] = F_d[kp + 7*size];

		// c_9
		xp = x; 
		if(y < NY1) yp = y;
		else yp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 9*size] = F_d[kp + 9*size];

		// c_12
		yp = y; 
		if(x > N16) xp = x;
		else xp = N16+NX1;
		kp = GID(zp, yp, xp);
		F_d[k + 12*size] = F_d[kp + 12*size];

		// c_13
		yp = y; 
		if(x < N16+NX) xp = x;
		else xp = N16-1;
		kp = GID(zp, yp, xp);
		F_d[k + 13*size] = F_d[kp + 13*size];


		//upper
		z = NZ1;
		k = GID(z, y, x);
		// c_6
		zp = 0;
		yp = y; xp = x;
		kp = GID(zp, yp, xp);
		F_d[k + 6*size] = F_d[kp + 6*size];

		// c_8
		xp = x; 
		if(y < NY1) yp = y;
		else yp = 0;
		kp = GID(zp, yp, xp);
		F_d[k + 8*size] = F_d[kp + 8*size];

		// c_10
		xp = x; 
		if(y > 1) yp = y;
		else yp = NY1+1;
		kp = GID(zp, yp, xp);
		F_d[k + 10*size] = F_d[kp + 10*size];

		// c_11
		yp = y; 
		if(x < N16+NX) xp = x;
		else xp = N16-1;
		kp = GID(zp, yp, xp);
		F_d[k + 11*size] = F_d[kp + 11*size];

		// c_14
		yp = y; 
		if(x > N16) xp = x;
		else xp = N16+NX1;
		kp = GID(zp, yp, xp);
		F_d[k + 14*size] = F_d[kp + 14*size];
	}

}
