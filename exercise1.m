close all 
clear

Gamma = 3;
xc = 0.75;
yc = 0.5;

nx = 51; 
ny = 41;
xmin = -2.5;
ymin = -2;
xmax = 2.5;
ymax = 2;

xm = zeros(nx,ny);
ym = zeros(nx,ny);
psi = zeros(nx,ny);

for i=1:nx
    for j=1:ny
    
    xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
    ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
    psi(i,j) = psipv(xc,yc,Gamma,xm(i,j),ym(i,j));
    
    end 
end 

c = -0.4 : 0.2 : 1.2;
contour(xm,ym,psi,c)