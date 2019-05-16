%Week 2 Exercise 1

close all
clear

ReL = 2500;
x = linspace(0,1,101);
ue = ones(1,length(x));

Mthick = zeros(1,length(x));

Mthick2 = 0;

for xi = 2:length(x)
    
    Mthick2 = Mthick2+0.45*(ue(xi)^-6)*ueintbit(x(xi-1),ue(xi-1),x(xi),ue(xi))/ReL; 
    Mthick(xi) = sqrt(Mthick2);
    
end

MthickB = 0.664*x.^0.5/(ReL^0.5);

figure
plot(x,Mthick,x,MthickB)
legend('Mthick', 'MthickB')



