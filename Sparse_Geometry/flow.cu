#include "lb.h"
#include "common.h"

//used for MRT-LBE model
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
__global__ void Evol_flow(double rgama, double *sf_d, double dt, double Fx, double Fy, double Fz, int N, int size, int *node_d, double *f_d, double *F_d)
{
	double P, U, V, W;
	int tx;
	int bx;
	int k;	
	
	double mf[Q];
	double f0, f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18;
	double F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18;
	int Id_0, Id_1, Id_2, Id_3, Id_4, Id_5, Id_6, Id_7, Id_8, Id_9, Id_10, Id_11, Id_12, Id_13, Id_14, Id_15, Id_16, Id_17, Id_18;

	tx = threadIdx.x;
    bx = blockIdx.x;
	k = bx*BX+tx;


	if(k < N)
	{

		f0  = f_d[k+0*size]; 
		f1  = f_d[k+1*size]; 
		f2  = f_d[k+2*size]; 
		f3  = f_d[k+3*size]; 
		f4  = f_d[k+4*size]; 
		f5  = f_d[k+5*size]; 
		f6  = f_d[k+6*size]; 
		f7  = f_d[k+7*size]; 
		f8  = f_d[k+8*size]; 
		f9  = f_d[k+9*size];

		f10 = f_d[k+10*size]; 
		f11 = f_d[k+11*size]; 
		f12 = f_d[k+12*size]; 
		f13 = f_d[k+13*size]; 
		f14 = f_d[k+14*size]; 
		f15 = f_d[k+15*size]; 
		f16 = f_d[k+16*size]; 
		f17 = f_d[k+17*size]; 
		f18 = f_d[k+18*size];

		Id_0 = node_d[k + 0*size];
		Id_1 = node_d[k + 1*size];
		Id_2 = node_d[k + 2*size];
		Id_3 = node_d[k + 3*size];
		Id_4 = node_d[k + 4*size];
		Id_5 = node_d[k + 5*size];
		Id_6 = node_d[k + 6*size];
		Id_7 = node_d[k + 7*size];
		Id_8 = node_d[k + 8*size];
		Id_9 = node_d[k + 9*size];

		Id_10 = node_d[k + 10*size];
		Id_11 = node_d[k + 11*size];
		Id_12 = node_d[k + 12*size];
		Id_13 = node_d[k + 13*size];
		Id_14 = node_d[k + 14*size];
		Id_15 = node_d[k + 15*size];
		Id_16 = node_d[k + 16*size];
		Id_17 = node_d[k + 17*size];
		Id_18 = node_d[k + 18*size];
		
		
		//f-mf///////////////////////////
		mf[0]  = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18; 
 		mf[1]  = -30.0*f0 - 11.0*f1 - 11.0*f2 - 11.0*f3 - 11.0*f4 - 11.0*f5 -11.0*f6 + 8.0*f7 + 8.0*f8 + 8.0*f9 
				+ 8.0*f10 + 8.0*f11 + 8.0*f12 + 8.0*f13 + 8.0*f14 + 8.0*f15 + 8.0*f16 + 8.0*f17 + 8.0*f18;
        mf[2]  = 12.0*f0 - 4.0*f1 - 4.0*f2 - 4.0*f3 - 4.0*f4 - 4.0*f5 - 4.0*f6 + f7 + f8 + f9 
				+ f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
        mf[3]  = f1 - f2 - f11 + f12 - f13 + f14 - f15 + f16 - f17 + f18;
        mf[4]  = -4.0*f1 + 4.0*f2 - f11 + f12 - f13 + f14 - f15 + f16 - f17 + f18;
        mf[5]  = f3 - f4 + f7 - f8 - f9 + f10 + f15 - f16 - f17 + f18;
        mf[6]  = -4.0*f3 + 4.0*f4 + f7 - f8 - f9 + f10 + f15 - f16 - f17 + f18;
        mf[7]  = f5 - f6 + f7 - f8 + f9 + f12 - f11 - f10 + f13 - f14;
        mf[8]  = f12 - f11 - f10 + f13 - f14 - 4.0*f5 + 4.0*f6 + f7 - f8 + f9;
        mf[9]  = 2.0*f1 - 2.0*f10 + f11 + f12 + f13 + f14 + f15 + f16 
				+ f17 + f18 + 2.0*f2 - f3 - f4 - f5 - f6 - 2.0*f7 - 2.0*f8 - 2.0*f9;

        mf[10] = f11 - 2.0*f10 - 4.0*f1 + f12 + f13 + f14 + f15 
				+ f16 + f17 + f18 - 4.0*f2 + 2.0*f3 + 2.0*f4 + 2.0*f5 + 2.0*f6 - 2.0*f7 - 2.0*f8 - 2.0*f9;
        mf[11] = f15 - f12 - f13 - f14 - f11 + f16 + f17 + f18 + f3 + f4 - f5 - f6;
        mf[12] = f15 - f12 - f13 - f14 - f11 + f16 + f17 + f18 - 2.0*f3 - 2.0*f4 + 2.0*f5 + 2.0*f6;
        mf[13] = f17 - f16 - f15 + f18;
        mf[14] = f7 - f10 + f8 - f9;
        mf[15] = f11 + f12 - f13 - f14;
        mf[16] = f11 - f12 + f13 - f14 - f15 + f16 - f17 + f18;
        mf[17] = f10 - f15 + f16 + f17 - f18 + f7 - f8 - f9;
        mf[18] = f10 - f11 + f12 + f13 - f14 - f7 + f8 - f9;

		//macroscopic variables/////////////////////////////////////////////////////////////////////////////
		U   = f1 + f12 + f14 + f16 + f18 - f2 - f11 - f13 - f15 - f17 + 0.5*dt*Fx;
		V   = f3 + f7  + f10 + f15 + f18 - f4 - f8  - f9  - f16 - f17 + 0.5*dt*Fy;
		W   = f5 + f7  + f9  + f12 + f13 - f6 - f8  - f10 - f11 - f14 + 0.5*dt*Fz;
		P = (f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18);
		
		//collision//-------------------------------------------------------------------------------------------------------------------------------
		mf[0]  = (mf[0] - sf_d[0]*(mf[0] - MEQ_0(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[0])*F_0(U,V,W,Fx,Fy,Fz,rgama);
		mf[1]  = (mf[1] - sf_d[1]*(mf[1] - MEQ_1(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[1])*F_1(U,V,W,Fx,Fy,Fz,rgama);
		mf[2]  = (mf[2] - sf_d[2]*(mf[2] - MEQ_2(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[2])*F_2(U,V,W,Fx,Fy,Fz,rgama);
		mf[3]  = (mf[3] - sf_d[3]*(mf[3] - MEQ_3(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[3])*F_3(U,V,W,Fx,Fy,Fz,rgama);
		mf[4]  = (mf[4] - sf_d[4]*(mf[4] - MEQ_4(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[4])*F_4(U,V,W,Fx,Fy,Fz,rgama);
		mf[5]  = (mf[5] - sf_d[5]*(mf[5] - MEQ_5(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[5])*F_5(U,V,W,Fx,Fy,Fz,rgama);
		mf[6]  = (mf[6] - sf_d[6]*(mf[6] - MEQ_6(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[6])*F_6(U,V,W,Fx,Fy,Fz,rgama);
		mf[7]  = (mf[7] - sf_d[7]*(mf[7] - MEQ_7(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[7])*F_7(U,V,W,Fx,Fy,Fz,rgama);
		mf[8]  = (mf[8] - sf_d[8]*(mf[8] - MEQ_8(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[8])*F_8(U,V,W,Fx,Fy,Fz,rgama);
		mf[9]  = (mf[9] - sf_d[9]*(mf[9] - MEQ_9(P, U, V, W, rgama)))+dt*(1.0 - 0.5*sf_d[9])*F_9(U,V,W,Fx,Fy,Fz,rgama);

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
		F0  = r0*mf[0] - r1*mf[1] + r2*mf[2];
        F1  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[3] - r5*mf[4] + r6*mf[9] - r6*mf[10];
        F2  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[3] + r5*mf[4] + r6*mf[9] - r6*mf[10]; 
        F3  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[5] - r5*mf[6] - r7*mf[9] + r7*mf[10] + r11*mf[11] - r11*mf[12]; 
        F4  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[5] + r5*mf[6] - r7*mf[9] + r7*mf[10] + r11*mf[11] - r11*mf[12];
        F5  = r0*mf[0] - r3*mf[1] - r4*mf[2] + r5*mf[7] - r5*mf[8] - r7*mf[9] + r7*mf[10] - r11*mf[11] + r11*mf[12];
        F6  = r0*mf[0] - r3*mf[1] - r4*mf[2] - r5*mf[7] + r5*mf[8] - r7*mf[9] + r7*mf[10] - r11*mf[11] + r11*mf[12];
        F7  = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[5] + r10*mf[6] + r5*mf[7] + r10*mf[8] - r6*mf[9]  - r7*mf[10] + r12*mf[14] + r15*mf[17] - r15*mf[18];
        F8  = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[5] - r10*mf[6] - r5*mf[7] - r10*mf[8] - r6*mf[9]  - r7*mf[10] + r12*mf[14] - r15*mf[17] + r15*mf[18];
        F9  = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[5] - r10*mf[6] + r5*mf[7] + r10*mf[8] - r6*mf[9]  - r7*mf[10] - r12*mf[14] - r15*mf[17] - r15*mf[18];

        F10 = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[5] + r10*mf[6] - r5*mf[7] - r10*mf[8] - r6*mf[9]  - r7*mf[10] - r12*mf[14] + r15*mf[17] + r15*mf[18];
        F11 = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] - r5*mf[7] - r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] + r12*mf[15] + r15*mf[16] - r15*mf[18];
        F12 = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] + r5*mf[7] + r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] + r12*mf[15] - r15*mf[16] + r15*mf[18];
        F13 = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] + r5*mf[7] + r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] - r12*mf[15] + r15*mf[16] + r15*mf[18]; 
        F14 = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] - r5*mf[7] - r10*mf[8] + r7*mf[9] + r13*mf[10] - r11*mf[11] - r14*mf[12] - r12*mf[15] - r15*mf[16] - r15*mf[18]; 
        F15 = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] + r5*mf[5] + r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] - r12*mf[13] - r15*mf[16] - r15*mf[17]; 
        F16 = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] - r5*mf[5] - r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] - r12*mf[13] + r15*mf[16] + r15*mf[17]; 
        F17 = r0*mf[0] + r8*mf[1] + r9*mf[2] - r5*mf[3] - r10*mf[4] - r5*mf[5] - r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] + r12*mf[13] - r15*mf[16] + r15*mf[17];
        F18 = r0*mf[0] + r8*mf[1] + r9*mf[2] + r5*mf[3] + r10*mf[4] + r5*mf[5] + r10*mf[6] + r7*mf[9] + r13*mf[10] + r11*mf[11] + r14*mf[12] + r12*mf[13] + r15*mf[16] - r15*mf[17];
		
		//  streaming  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		F_d[Id_0] = F0;
		F_d[Id_1] = F1;
		F_d[Id_2] = F2;
		F_d[Id_3] = F3;
		F_d[Id_4] = F4;
		F_d[Id_5] = F5;
		F_d[Id_6] = F6;
		F_d[Id_7] = F7;
		F_d[Id_8] = F8;
		F_d[Id_9] = F9;

		F_d[Id_10] = F10;
		F_d[Id_11] = F11;
		F_d[Id_12] = F12;
		F_d[Id_13] = F13;
		F_d[Id_14] = F14;
		F_d[Id_15] = F15;
		F_d[Id_16] = F16;
		F_d[Id_17] = F17;
		F_d[Id_18] = F18;
	
	}
}
