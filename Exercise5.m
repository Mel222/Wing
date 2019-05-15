close all
clear

np=100;

theta=(0:np)*2*pi/np;

xs=cos(theta);    %array of x coordinates

ys=sin(theta);  %array of y coordinates

alpha=0;
A0 = build_lhs(xs,ys);
b0 = build_rhs(xs,ys,alpha);
gamma0 = A0\b0;

figure(10)
plot(theta/pi, gamma0)
axis([0 2 -2.5 2.5])

alpha=0.1;
A1 = build_lhs(xs,ys);
b1 = build_rhs(xs,ys,alpha);
gamma1 = A1\b1;

figure(11)
plot(theta/pi, gamma1)
axis([0 2 -2.5 2.5])

Circ0 = sum(gamma0)*theta(2)
Circ1 = sum(gamma1)*theta(2)

CircTheo = -4*pi*sin(alpha)

