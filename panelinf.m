
function [infa, infb]= panelinf(xa,ya,xb,yb,x,y)

del=sqrt((yb-ya)^2+(xb-xa)^2);

% unit tangential vector to the panel
tx=1/sqrt((xb-xa)^2+(yb-ya)^2)*(xb-xa);
ty=1/sqrt((xb-xa)^2+(yb-ya)^2)*(yb-ya);
% unit normal vector to the panel
nx=1/sqrt((xb-xa)^2+(yb-ya)^2)*(ya-yb);
ny=1/sqrt((xb-xa)^2+(yb-ya)^2)*(xb-xa);

rx=(x-xa);
ry=(y-ya);

X=rx*tx+ry*ty;
Y=rx*nx+ry*ny;

[infa,infb]=refpaninf(del,X,Y);
end

