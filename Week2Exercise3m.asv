%Week 2 Exercise 3

close all
clear

% ReL = [1e05, 1e04, 1e03]; %% three dif 
ReL = linspace(1.78e06, 1.8e06, 1000); 
x = linspace(0,1,101);

uegrad = -0.5;
ue = zeros(1,length(x));
ue(1) = 1;

for i = 2:length(x)
    ue(i) = ue(i-1) + uegrad*x(2);
end

figure 
plot(x,ue)

int = zeros(1,length(ReL));
ils = zeros(1,length(ReL));

Supplant = true;

for r = 1:length(ReL)
        
    integral = 0; 
    Mthick = 0;

    n = 101;
    laminar = true;
    i = 1;

    while laminar && i < n

        Rethet = ReL(r)*ue(i)*Mthick;
        m = -ReL(r)*(Mthick^2)*uegrad;

        H = thwaites_lookup(m);
        He = laminar_He(H);

        if log(Rethet) >= 18.4*He - 21.74
            laminar = false;
%             disp([x(i) Rethet/1000])
            int(r) = x(i);
        elseif m >= 0.09
            laminar = false;
%             disp([x(i) Rethet/1000])
            ils(r) = x(i);
        end

        integral = integral + ueintbit(x(i),ue(i),x(i+1),ue(i+1)); 
        Mthick = sqrt(0.45*(ue(i+1)^-6)*integral/ReL(r));

        i = i+1;
    end
    
    if ils(r)==0 && Supplant
        Supplant = false; 
        disp(['Transition supplants seperation ' num2str(ReL(r))])
    end
        
    if int(r)~= 0
%         disp(['Natural transition at ' num2str(int(r)) ' with Rethet ' num2str(Rethet)]) 
    elseif ils(r)~= 0
%         disp(['Laminar separation at ' num2str(ils(r)) ' with Rethet ' num2str(Rethet)])
    end    
end 

% 1791352.7914
% 1791351.3514

