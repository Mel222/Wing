close all
clear

%Definition of the cylinder panels

np=100;

theta=(0:np)*2*pi/np;

xs=cos(theta);    %array of x coordinates

ys=sin(theta);  %array of y coordinates

gamma=-2*sin(theta); %array of the vortex sheet strength



nx = 51;

ny = 41;

xmin = -2.5;

ymin = -2;

xmax = 2.5;

ymax = 2;



xm = zeros(nx,ny);

ym = zeros(nx,ny);

psiU=zeros(nx,ny);



for i=1:nx

    for j=1:ny

        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);

        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);

        psiU(i,j)=ym(i,j);  %setting the free stream contribution to the stream function

    end    

end        

psi = zeros(nx,ny);   

for i=1:nx

    for j=1:ny
        
        psiTOT=0;
        
        for p =1:np
   
            [infaP, infbP] = panelinf(xs(p),ys(p),xs(p+1),ys(p+1),xm(i,j),ym(i,j));
            psiP = gamma(p)*infaP + gamma(p+1)*infbP;
            psiTOT = psiTOT + psiP;
            
        end
        
        psi(i,j)= psiTOT;

    end 

end 

psiFinal = psi+psiU;

figure(9)

c = -1.75 : 0.25 : 1.75;

contour(xm,ym,psiFinal,c)

hold on
plot(xs,ys)
hold off