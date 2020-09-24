clc; clear;

Datafile='flow.dat';

formatstr='%d %d %d %f %f %f %f %d';

NX=31;

NY=NX;

NZ=NX;

file=fopen(Datafile);
Data=textscan(file,formatstr,'Headerlines',3);                         
ux=Data{1,4};
uy=Data{1,5};
uz=Data{1,6};
P=Data{1,7};
phi=Data{1,8};
fclose(file);

u=zeros(NZ,NY,NX);
v=zeros(NZ,NY,NX);
w=zeros(NZ,NY,NX);
p=zeros(NZ,NY,NX);
% flag=ones(NZ,NY,NX);

for z=1:NZ;
    for y=1:NY;
        for x=1:NX;
            u(z,y,x)=ux((z-1)*NX*NY+(y-1)*NX+x);
            v(z,y,x)=uy((z-1)*NX*NY+(y-1)*NX+x);
            w(z,y,x)=uz((z-1)*NX*NY+(y-1)*NX+x);
            p(z,y,x)=P((z-1)*NX*NY+(y-1)*NX+x);
            flag(z,y,x)=phi((z-1)*NX*NY+(y-1)*NX+x);
        end
    end
end

sum = 0;

for z=1:NZ;
    for y=1:NY;
        for x=1:NX;
            sum = sum + u(z,y,x);
        end
    end
end

au = sum/NX/NY/NZ;

pin = p(:,:,1);

pout = p(:,:,NX);

fg = flag(:,:,1);

sumpin = 0.0;

sumpout = 0.0;

N = 0;

for z = 1:NZ;
    for y = 1:NY;
        if fg(z,y) == 1
            sumpin = sumpin + pin(z,y);
            sumpout= sumpout + pout(z,y);
            N = N + 1;
        end
    end
end

dp = (sumpin - sumpout)/N;

nu = 0.1/(NX-1);

D = 0.7;

Re = 0.1;

U = Re*nu/D;

Ka = 0.025448216331798;

G = U*nu/Ka;

K = nu*au/G

% for z=1:NZ;
%     for y=1:NY;
%         x = 1;
%         ui(z,y)=ux((z-1)*NX*NY+(y-1)*NX+x);
%         x = NX;
%         uo(z,y)=ux((z-1)*NX*NY+(y-1)*NX+x);
%     end
% end

% save uc.dat uc -ascii;

% figure;
% mesh(ui);
% 
% figure;
% mesh(uo);




            
            

