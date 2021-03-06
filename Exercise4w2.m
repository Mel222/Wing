clear 
close all

global Re ue0 duedx
Re=1e7;
ue0=1;  %free stream velocity
duedx=0;

%Definition of the initial value of theta and deltaE
x0 = 0.01;
thick0(1)=0.037*x0*(Re*x0)^(-1/5);
thick0(2)=1.80*thick0(1);

%the differential equation is solved

[delx, thickhist]=ode45(@thickdash,[0 0.99], thick0);

x=x0+delx

theta7=(0.037*x.*(Re*x).^(-1/5))  % .* operation done element by element of the vector (elemnetwise op??)
theta9=(0.023*x.*(Re*x).^(-1/6))

plot(x,thickhist(:,1),x,theta7,x,theta9)

legend( 'theta', 'theta7', 'theta9')



