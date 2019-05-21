%Week 2 Exercise 6

close all
clear

global Re ue0 duedx

ReL = 1e5; 
x = linspace(0,1,101);

uegrad= linspace(-0.38125,-0.37773,5);
for u = 1:length(uegrad)
duedx = uegrad(u);
ue = zeros(1,length(x));
ue(1) = 1;
for i = 2:length(x)
    ue(i) = ue(i-1) + uegrad(u)*x(2);
end

He=zeros(1,length(x));
He(1)=1.57258; %an arbitrary value is assined to He at point x=0
Theta=zeros(1,length(x));

% figure 
% plot(x,ue)

% int = zeros(1,length(ReL)); %location of natural transition
% ils = zeros(1,length(ReL)); %location of natural separation
% itr = zeros(1,length(ReL)); %location of turbulent reattachement
% its = zeros(1,length(ReL)); %location of turbulent separation

int = 0; %location of natural transition
ils = 0;%location of natural separation
itr = 0; %location of turbulent reattachement
its = 0; %location of turbulent separation

%laminar loop


        
integral = 0; 
Mthick = 0;

n = 101;
laminar = true;
i = 1;

while laminar && i < n
    i = i + 1;
    
    integral = integral + ueintbit(x(i-1),ue(i-1),x(i),ue(i)); 
    Mthick = sqrt(0.45*(ue(i)^-6)*integral/ReL);
    Theta(i)=Mthick;
    
    Rethet = ReL*ue(i)*Mthick;
    m = -ReL*(Mthick^2)*uegrad(u);

    H = thwaites_lookup(m);
    HE = laminar_He(H);

    He(i)=HE;    
    if log(Rethet) >= 18.4*HE - 21.74
        laminar = false;
%             disp([x(i) Rethet/1000])
        int = i;
    elseif m >= 0.09
        laminar = false;  %if laminar separation is detected
%             disp([x(i) Rethet/1000])
        ils = i;
        He(i)=1.51509;   %when laminar separation is detected, He is set to its laminar separation value
    end
        
end
        
    if int~= 0
        disp(['Natural transition at ' num2str(int) ' with Rethet ' num2str(Rethet) ' with uegrad ' num2str(uegrad(u))]) 
    elseif ils~= 0
        disp(['Laminar separation at ' num2str(ils) ' with Rethet ' num2str(Rethet) ' with uegrad ' num2str(uegrad(u))])
    end   

%step (iii)
    if int~=0  %ie natural transition occured
       deltaE=He(int)*Theta(int); %deltaE is computed where transition occured
    elseif ils~=0 %ie laminar separtion occured
        deltaE=He(ils)*Theta(ils);
    else
        deltaE=He(i-1)*Theta(i-1);
    end



%while loop for the turbulent boundary layer

thick0(1) = Theta(i);
thick0(2) = deltaE;

while i<length(x) && its==0
    i = i+1;
    Re = ReL;
    ue0 = ue(i);
    
    [delx, thickhist] = ode45(@thickdash,[0,x(i) - x(i-1)],thick0);
        
     thick0(1) = thickhist(end,1);   
     thick0(2) = thickhist(end,2); 
     HE = thick0(2)/thick0(1);
     Rethet = ReL*ue(i)*thick0(1);
     
     Theta(i) = thick0(1);
     He(i) = HE;
     
     if HE <= 1.46
         its = i;
         disp(['Turbulent separation at ' num2str(its) ' with Rethet ' num2str(Rethet) ' with uegrad ' num2str(uegrad(u))])
  
     elseif ils ~= 0 && HE >= 1.58 && itr == 0 
         itr = i;
         disp(['Turbulent reattachement at ' num2str(itr) ' with Rethet ' num2str(Rethet) ' with uegrad ' num2str(uegrad(u))])
     end
     
end

if its ~= 0
    H = 2.803;
    for i = its:(length(x)-1)
        Theta(i+1) = Theta(i)*(ue(i)/ue(i+1))^(H+2);
    end
end
if its == length(x)
    disp(['Seperates at the end ' num2str(duedx)])
end
end

% figure(2)
% plot(x,Theta)
% figure(3)
% plot(x,He)

% Seperates at end for uegrad = -0.38 
