close all
clear
% Plot of the exact influence coefficients
nx = 51;
ny = 41;
xmin = 0;
ymin = 0;
xmax = 5;
ymax = 4;
xa=4.1;
xb=2.2;
ya=1.3;
yb=2.9;


xm = zeros(nx,ny);
ym = zeros(nx,ny);
infa3 = zeros(nx,ny);
infb3 = zeros(nx,ny);

for i=1:nx
    for j=1:ny
    
    xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
    ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
    [infa3(i,j), infb3(i,j)]=panelinf(xa,ya,xb,yb,xm(i,j),ym(i,j))
    
    end 
end 

figure(5)
c = -0.15 : 0.05 : 0.15;
contour(xm,ym,infa3,c)
figure(6)
contour(xm,ym,infb3,c)

%Plot of the approximate influence coefficients


nx = 51; 
ny = 41;
nv = 100;
xmin = 0;
ymin = 0;
xmax = 5;
ymax = 4;
xa=4.1;
xb=2.2;
ya=1.3;
yb=2.9;
del=sqrt((xb-xa)^2+(yb-ya)^2)


xm = zeros(nx,ny);
ym = zeros(nx,ny);
infa3approx = zeros(nx,ny);
infb3approx=zeros(nx,ny);


for i=1:nx
    for j=1:ny
    xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
    ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);

    infa3TOT = 0;
    infb3TOT = 0;

       for v=1:nv

        xc = xa+(del/(2*nv)+(v-1)*del/nv)*1/sqrt((xb-xa)^2+(yb-ya)^2)*(xb-xa);
        yc = ya+(yb-ya)/(xb-xa)*(xc-xa);
        infa3N = psipv(xc,yc,(1-sqrt((xc-xa)^2+(yc-ya)^2)/del)*del/nv,xm(i,j),ym(i,j));
        %(1-sqrt((xc-xa)^2+(yc-ya)^2)/del)*del/nv strength of the vortex associated with
        %the influence coefficient fa
        infb3N= psipv(xc,yc,sqrt((xc-xa)^2+(yc-ya)^2)/del*del/nv,xm(i,j),ym(i,j));
        % sqrt((xc-xa)^2+(yc-ya)^2)/del*del/nv is the strength of the vortex
        %associated with the influence coefficient fb
        infa3TOT= infa3TOT + infa3N;
        infb3TOT= infb3TOT + infb3N;
       end
        infa3approx(i,j) = infa3TOT;
        infb3approx(i,j) = infb3TOT;
    end 
end 


figure(7)
c = -0.15 : 0.05 : 0.15;
contour(xm,ym,infa3approx,c)
figure(8)
contour(xm,ym,infb3approx,c)
