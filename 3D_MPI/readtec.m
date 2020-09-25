clear all;

Datafile='TPL6000_0.dat';

formatstr='%d %d %d %f %f %f %f %d';
NX=25;NY=50;NZ=50;
file=fopen(Datafile);
Data=textscan(file,formatstr,'Headerlines',3);                         
ux=Data{1,4};
uy=Data{1,5};
uz=Data{1,6};
rho=Data{1,7};
phi=Data{1,8};
fclose(file);

u=zeros(NZ,NY,NX);
v=zeros(NZ,NY,NX);
w=zeros(NZ,NY,NX);
p=zeros(NZ,NY,NX);
flag=zeros(NZ,NY,NX);

for z=1:NZ;
    for y=1:NY;
        for x=1:NX;
            u(z,y,x)=ux((z-1)*NX*NY+(y-1)*NX+x);
            v(z,y,x)=uy((z-1)*NX*NY+(y-1)*NX+x);
            w(z,y,x)=uz((z-1)*NX*NY+(y-1)*NX+x);
            p(z,y,x)=rho((z-1)*NX*NY+(y-1)*NX+x)/3;
            flag(z,y,x)=phi((z-1)*NX*NY+(y-1)*NX+x);
        end
    end
end

for z=1:NZ;
    for y=1:NY;
        x = 1;
        ui(z,y)=ux((z-1)*NX*NY+(y-1)*NX+x);
        x = NX;
        uo(z,y)=ux((z-1)*NX*NY+(y-1)*NX+x);
    end
end

% save uc.dat uc -ascii;

figure;
mesh(ui);

figure;
mesh(uo);




            
            

