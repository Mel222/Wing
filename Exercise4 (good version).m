

%Definition of the cylinder panels

np=100
theta=(0:np)*2*pi/np
xs=cos(theta)    %array of x coordinates
ys=sin(theta)  %array of y coordinates
gamma=-2*sin(theta) %array of the vortex sheet strength

nx = 51;
ny = 41;
xmin = -2.5;
ymin = -2;
xmax = 2.5;
ymax = 2;

xm = zeros(nx,ny);
ym = zeros(nx,ny);
psi=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psi(i,j)=ym(i,j)  %setting the free stream contribution to the stream function
    end    
end        
   
for i=1:nx
    for j=1:ny
    
    
   [infa3(i,j), infb3(i,j)]=panelinf(xs(i),ys(i),xs(i+1),ys(i+1),xm(i,j),ym(i,j))
    
    end 
end 