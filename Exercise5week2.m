clear 
close all

global Re ue0 duedx
Re=1e7;
ue0=1;  %free stream velocity
duedx=-0.6


%Definition of the initial value of theta and deltaE
x0=0.01;
thick0(1)=0.037*x0*(Re*x0)^(-1/5);
thick0(2)=1.80*thick0(1);

%the differential equation is solved

[delx thickhist]=ode45(@thickdash,[0 0.99], thick0);
size(thickhist)

x=x0+delx;

HE=thickhist(:,2)./thickhist(:,1);

seplocation=0;
i=1;

 while i<length(HE) && seplocation==0
    if HE(i)<=1.46
        seplocation=x(i)
    end
    i=i+1;
 end   
 


figure(1)
plot(x, HE)

figure(2)
plot(x, thickhist(:,1),x,thickhist(:,2))
legend('theta','deltaE')