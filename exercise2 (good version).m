close all
clear 
% Plot of the exact influence coefficients
nx = 51;
ny = 41;
xmin = -2.5;
ymin = -2;
xmax = 2.5;
ymax = 2;

del = 1.5;

xm = zeros(nx,ny);
ym = zeros(nx,ny);
infa = zeros(nx,ny);
infb = zeros(nx,ny);

for i=1:nx
    for j=1:ny
    
    xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
    ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
    [infa(i,j), infb(i,j)] = refpaninf(del,xm(i,j),ym(i,j));
    
    end 
end 

figure(1)
c = -0.15 : 0.05 : 0.15;
contour(xm,ym,infa,c)
figure(2)
contour(xm,ym,infb,c)


%Plot of the approximate influence coefficients


nx = 51; 
ny = 41;
nv = 100;
xmin = -2.5;
ymin = -2;
xmax = 2.5;
ymax = 2;


xm = zeros(nx,ny);
ym = zeros(nx,ny);
infaapprox = zeros(nx,ny);
infbapprox=zeros(nx,ny)


for i=1:nx
    for j=1:ny
    xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
    ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);

    infaTOT = 0;
    infbTOT = 0;

       for v=1:nv

        xc = del/(nv*2) + (v-1)*(del)/(nv-1);
        yc = 0;
        infaN = psipv(xc,yc,(1-xc/del)*del/nv,xm(i,j),ym(i,j));
        %(1-xc/del)*del/nv strength of the vortex associated with the influence coefficient fa
        infbN= psipv(xc,yc,xc/del*del/nv,xm(i,j),ym(i,j));
        %9xc/del)*del/nv strength of the vortex associated with the influence coefficient fb
        infaTOT= infaTOT + infaN;
        infbTOT= infbTOT + infbN;
       end
        infaapprox(i,j) = infaTOT;
        infbapprox(i,j) = infbTOT;
    end 
end 


figure(3)
c = -0.15 : 0.05 : 0.15;
contour(xm,ym,infaapprox,c)
figure(4)
contour(xm,ym,infbapprox,c)
