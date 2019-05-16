%Week 2 Exercise 1

close all
clear

ReL = 10e06; %% three dif 
x = linspace(0,1,101);

uegrad = -0.2; %% three dif 
ue = zeros(1,length(x));
ue(1) = 1;

for i = 2:length(ue)
    ue(i) = ue(i-1) + uegrad*x(2);
end

figure 
plot(x,ue)

%% Vector Method 

% Mthick = zeros(1,length(x));
% 
% Mthick2 = 0;
% 
% for xi = 2:length(x)
%     
%     Mthick2 = Mthick2 + 0.45*(ue(xi)^-6)*ueintbit(x(xi-1),ue(xi-1),x(xi),ue(xi))/ReL; 
%     Mthick(xi) = sqrt(Mthick2);
%     
% end
% 
% Rethet = ReL.*ue.*Mthick;
% m = -ReL.*(Mthick.^2)*uegrad;
% 
% H = arrayfun(@thwaites_lookup,m);
% He = arrayfun(@laminar_He,H);
% 
% n = 101;
% laminar = true;
% i = 1;
% while laminar && i < n
% %     Rethet(i)
%     if log(Rethet(i)) >= 18.4*He(i) - 21.74
%         laminar = false;
%         disp([x(i) Rethet(i)/1000])
%         
%     end
%     i = i+1;
% end 

%% Iterative Method 

integral = 0; 
Mthick = 0;

% for i = 1:(length(x)-1)
%     
% %     Mthick2 = Mthick2 + 0.45*(ue(i)^-6)*ueintbit(x(i),ue(i),x(i+1),ue(i+1))/ReL; 
% %     Mthick = sqrt(Mthick2);
%     
%     Rethet = ReL*ue(i)*Mthick;
%     m = -ReL*(Mthick^2)*uegrad;
% 
%     H = thwaites_lookup(m);
%     He = laminar_He(H);
%     
%     if log(Rethet) >= 18.4*He - 21.74
%         laminar = false;
% %         disp([x(i) Rethet/1000])
%     end
%         
%     integral = integral + ueintbit(x(i),ue(i),x(i+1),ue(i+1)); 
%     Mthick = sqrt(0.45*(ue(i+1)^-6)*integral/ReL);
%  
% end

n = 101;
laminar = true;
i = 1;

while laminar && i < n
    Rethet = ReL*ue(i)*Mthick;
    m = -ReL*(Mthick^2)*uegrad;

    H = thwaites_lookup(m);
    He = laminar_He(H);
    
    if log(Rethet) >= 18.4*He - 21.74
        laminar = false;
        disp([x(i) Rethet/1000])
    end
        
    integral = integral + ueintbit(x(i),ue(i),x(i+1),ue(i+1)); 
    Mthick = sqrt(0.45*(ue(i+1)^-6)*integral/ReL);
        
    i = i+1;
end
 





