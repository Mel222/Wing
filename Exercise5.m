close all
clear

np=100;

theta=(0:np)*2*pi/np;

xs=cos(theta);    %array of x coordinates

ys=sin(theta);  %array of y coordinates


nx = 51;

ny = 41;

xmin = -2.5;

ymin = -2;

xmax = 2.5;

ymax = 2;



xm = zeros(nx,ny);

ym = zeros(nx,ny);



for i=1:nx

    for j=1:ny

        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);

        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);


    end    

end

alpha=0;
A = build_lhs(xs,ys);
b = build_rhs(xs,ys,alpha);
gamma = A\b;

figure(10)
plot(gamma)

