#include "lb.h"
#include "common.h"


// save the results
void datasave()
{
	int z, y, x, k;
	int q;
	double ut, vt, wt, uvwt, pt;
	FILE *fp; 
	FILE *fp1, *fp2, *fp3, *fp4;

	if((fp=fopen("flow.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	if((fp1=fopen("u.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	if((fp2=fopen("v.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	if((fp3=fopen("w.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	if((fp4=fopen("p.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	fprintf(fp,"Title=\"LBM Driven Flows\"\n");
	fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\",\"U\",\"V\",\"W\",\"P\",\"flag\"\n");
	fprintf(fp,"ZONE I=%d,J=%d,K=%d,F=POINT\n",NX1,NY1,NZ1);
	for(z = 1; z <= NZ1; z++)
		for(y = 1; y <= NY1; y++)
			for(x = N16; x < N16+NX1; x++)
			{
				k = GID(z,y,x);
				if(flag[k] == 0)
				{
					ut = vt = wt = pt = 0.0;
					for(q = 0; q < Q; q++)
					{
						ut += e[q][0]*f[k+q*size];
						vt += e[q][1]*f[k+q*size];
						wt += e[q][2]*f[k+q*size];
						pt += f[k+q*size];
					}
					ut += 0.5*dt*Fx;
					vt += 0.5*dt*Fy;
					wt += 0.5*dt*Fz;
					uvwt = ut*ut + vt*vt + wt*wt;

				}
				else
				{
					ut = vt = wt = pt = 0.0;
				}
				fprintf(fp1, "%e ", ut);
				fprintf(fp2, "%e ", vt);
				fprintf(fp3, "%e ", wt);
				fprintf(fp4, "%e ", pt);
				fprintf(fp,"%d %d %d %e %e %e %e %d\n", x, y, z, ut, vt, wt, pt, flag[k]);
			}

	fclose(fp);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	x = N16;
	if((fp=fopen("uinlet.dat","w")) == NULL )  printf("FILE OPEN ERROR!\n");
//	if((fp1=fopen("pinlet.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	for(z = 1; z <= NZ1; z++)
	{
		for(y = 1; y <= NY1; y++)
		{
			ut = vt = wt = pt = 0.0;
			k = GID(z,y,x);
			if(flag[k] == 0)
			{
				for(q = 0; q < Q; q++)
				{
					ut += e[q][0]*f[k+q*size];
					vt += e[q][1]*f[k+q*size];
					wt += e[q][2]*f[k+q*size];
					pt += f[k+q*size];
				}
				ut += 0.5*dt*Fx;
				vt += 0.5*dt*Fy;
				wt += 0.5*dt*Fz;
				uvwt = ut*ut + vt*vt + wt*wt;
//				pt = 0.5*(pt - 0.5*uvwt);	
			}
			else
			{
				ut = vt = wt = pt = 0.0;
			}
			fprintf(fp, "%e ", ut);
//			fprintf(fp1, "%e ", pt);
		}
		fprintf(fp, "\n");
//		fprintf(fp1, "\n");
	}
	fclose(fp);
//	fclose(fp1);
/*
	x = N16+NX1/2;
	if((fp=fopen("uxc.dat","w")) == NULL )  printf("FILE OPEN ERROR!\n");
	if((fp1=fopen("pxc.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	for(z = 1; z <= NZ1; z++)
	{
		for(y = 1; y <= NY1; y++)
		{
			ut = vt = wt = pt = 0.0;
			k = GID(z,y,x);
			if(flag[k] == 1)
			{
				for(q = 1; q < Q; q++)
				{
					ut += e[q][0]*f[k+q*size];
					vt += e[q][1]*f[k+q*size];
					wt += e[q][2]*f[k+q*size];
					pt += f[k+q*size];
				}
				uvwt = ut*ut + vt*vt + wt*wt;
				pt = 0.5*(pt - 0.5*uvwt);	
			}
			else
			{
				ut = vt = wt = pt = 0.0;
			}
			fprintf(fp, "%e ", ut);
			fprintf(fp1, "%e ", pt);
		}
		fprintf(fp, "\n");
		fprintf(fp1, "\n");
	}
	fclose(fp);
	fclose(fp1);
*/
	x = N16+NX;
	if((fp=fopen("uoutlet.dat","w")) == NULL )  printf("FILE OPEN ERROR!\n");
//	if((fp1=fopen("poutlet.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	for(z = 1; z <= NZ1; z++)
	{
		for(y = 1; y <= NY1; y++)
		{
			ut = vt = wt = pt = 0.0;
			k = GID(z,y,x);
			if(flag[k] == 0)
			{
				for(q = 0; q < Q; q++)
				{
					ut += e[q][0]*f[k+q*size];
					vt += e[q][1]*f[k+q*size];
					wt += e[q][2]*f[k+q*size];
					pt += f[k+q*size];
				}
				ut += 0.5*dt*Fx;
				vt += 0.5*dt*Fy;
				wt += 0.5*dt*Fz;
				uvwt = ut*ut + vt*vt + wt*wt;
//				pt = 0.5*(pt - 0.5*uvwt);	
			}
			else
			{
				ut = vt = wt = pt = 0.0;
			}
			fprintf(fp, "%e ", ut);
//			fprintf(fp1, "%e ", pt);
		}
		fprintf(fp, "\n");
//		fprintf(fp1, "\n");
	}
	fclose(fp);
//	fclose(fp1);


	z = NZ2/2;
	if((fp=fopen("ucenter.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
//	if((fp1=fopen("pcenter.dat","w")) == NULL ) printf("FILE OPEN ERROR!\n");
	for(y = 1; y <= NY1; y++)
	{
		for(x = N16; x < N16+NX1; x++)
		{
			k = GID(z,y,x);
			if(flag[k] == 0)
			{
				ut = vt = wt = pt = 0.0;
				for(q = 0; q < Q; q++)
				{
					ut += e[q][0]*f[k+q*size];
					vt += e[q][1]*f[k+q*size];
					wt += e[q][2]*f[k+q*size];
					pt += f[k+q*size];
				}
				ut += 0.5*dt*Fx;
				vt += 0.5*dt*Fy;
				wt += 0.5*dt*Fz;
				uvwt = ut*ut + vt*vt + wt*wt;
//				pt = 0.5*(pt - 0.5*uvwt);	
			}
			else
			{
				ut = vt = wt = pt = 0.0;
			}
			fprintf(fp, "%e ", ut);
//			fprintf(fp1, "%e ", pt);
		}
		fprintf(fp, "\n");
//		fprintf(fp1, "\n");
	}
	fclose(fp);
//	fclose(fp1);

}
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
