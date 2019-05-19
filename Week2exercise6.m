%Week 2 Exercise 3

close all
clear

% ReL = [1e05, 1e04, 1e03]; %% three dif 
ReL = linspace(1.78e06, 1.8e06, 1000); 
x = linspace(0,1,101);

uegrad = -0.5;
ue = zeros(1,length(x));
ue(1) = 1;

He=zeros(1,lentgh(x));
He(0)=1.57258; %an arbitrary value is assined to He at point x=0
Theta=zeros(1,length(x))

for i = 2:length(x)
    ue(i) = ue(i-1) + uegrad*x(2);
end

figure 
plot(x,ue)

int = zeros(1,length(ReL)); %location of natural transition
ils = zeros(1,length(ReL)); %location of natural separation
itr = zeros(1,length(ReL)); %location of turbulent reattachement
its = zeros(1,length(ReL)); %location of turbulent separation

Supplant = true;

%laminar loop

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
        He(i)=He    
        if log(Rethet) >= 18.4*He - 21.74
            laminar = false;
%             disp([x(i) Rethet/1000])
            int(r) = i;
        elseif m >= 0.09
            laminar = false;  %if laminar separation is detected
%             disp([x(i) Rethet/1000])
            ils(r) = i;
            He(i)=1.51509;   %when laminar separation is detected, He is set to its laminar separation value
        end

        integral = integral + ueintbit(x(i),ue(i),x(i+1),ue(i+1)); 
        Mthick = sqrt(0.45*(ue(i+1)^-6)*integral/ReL(r));
        Theta(i)=Mthick;
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

%step (iii)
for r=1:3
    if int(r)~=0;  %ie natural transition occured
       deltaE=He(int(r))*Theta(int(r)); %deltaE is computed where transition occured
    elseif ils(r)~=0; %ie laminar separtion occured
        deltaE=He(int(r))*Theta(int(r));
    elseif
        deltaE=He(length(x))*Theta(length(x));
    end
end 


%while loop for the turbulent boundary layer

for r=1:length(ReL)
    
    if int(r)~=0
        i=int(r)
    elseif ils(r)~=0
        i=ils(r)
    elseif 
        i=length(x) % if there is no transition to turbulence
    end
    
    while i<length(x) && its(r)=0
        