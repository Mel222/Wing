function[int, ils, itr, its, delstar, theta] = bl_solv(x,cp)

global Re ue0 duedx
    
He=zeros(1,length(x));
He(1)=1.57258; %an arbitrary value is assined to He at point x=0
theta=zeros(1,length(x));
delstar=zeros(1,length(x)); 

int = 0; %location of natural transition
ils = 0;%location of natural separation
itr = 0; %location of turbulent reattachement
its = 0; %location of turbulent separation

%First Panel
        
integral = 0; 
n = length(x);
laminar = true;
i = 1;

ue(1) = (1-cp(1))^0.5;
uegrad1 = ue(1)/x(1);
duedx = uegrad1;
integral = integral + ueintbit(0,0,x(1),ue(1)); 
Mthick = sqrt(0.45*(ue(1)^-6)*integral/Re)
theta(1)=Mthick;
Rethet = Re*ue(1)*Mthick;
m = -Re*(Mthick^2)*uegrad1;
H = thwaites_lookup(m);
HE = laminar_He(H);
delstar(i)=H*theta(i);
He(1)=HE; 
    
if log(Rethet) >= 18.4*HE - 21.74
    laminar = false;
    int = i;
elseif m >= 0.09
    laminar = false;  %if laminar separation is detected
    ils = i;
    He(i)=1.51509;   %when laminar separation is detected, He is set to its laminar separation value
end
    
% Laminar loop

while laminar && i < n
    i = i + 1;
    
    ue(i) = (1-cp(i))^0.5;
    ue(i-1) = (1-cp(i-1))^0.5;
    
    uegrad = (ue(i)-ue(i-1))/(x(i)-x(i-1));
    duedx = uegrad;
    
    integral = integral + ueintbit(x(i-1),ue(i-1),x(i),ue(i))
    Mthick = sqrt(0.45*(ue(i)^-6)*integral/Re);
    theta(i)=Mthick;
    Mthick;
    
    Rethet = Re*ue(i)*Mthick;
    m = -Re*(Mthick^2)*uegrad;

    H = thwaites_lookup(m);
    HE = laminar_He(H);

    delstar(i)=H*theta(i);
    He(i)=HE;    
    if log(Rethet) >= 18.4*HE - 21.74
        laminar = false;
        int = i;
    elseif m >= 0.09
        laminar = false;  %if laminar separation is detected
        ils = i;
        He(i)=1.51509;   %when laminar separation is detected, He is set to its laminar separation value
    end
        
end     

deltaE=He(i)*theta(i);

thick0(1) = theta(i);
thick0(2) = deltaE;

%while loop for the turbulent boundary layer

while i<length(x) && its==0
    i = i+1;
    ue0 = ue(i-1);
    
    ue(i) = (1-cp(i))^0.5;
    ue(i-1) = (1-cp(i-1))^0.5;
    
    uegrad = (ue(i)-ue(i-1))/(x(i)-x(i-1));
    duedx = uegrad;
    
    [delx, thickhist] = ode45(@thickdash,[0,(x(i) - x(i-1))],thick0);
        
     thick0(1) = thickhist(end,1);   
     thick0(2) = thickhist(end,2); 
     HE = thick0(2)/thick0(1);
     Rethet = Re*ue(i)*thick0(1);
     
     theta(i) = thick0(1)
     He(i) = HE;
     
     H=(11*HE+15)/(48*HE-59);
     delstar(i)=H*theta(i);
     
     if HE < 1.46
         its = i;
  
     elseif ils ~= 0 && HE > 1.58 && itr == 0
         itr = i;
     end
     
end

if its ~= 0
    H = 2.803;
    for i = its:(length(x)-1)
        ue(i) = (1-cp(i))^0.5;
        ue(i+1) = (1-cp(i+1))^0.5;
        theta(i+1) = theta(i)*(ue(i)/ue(i+1))^(H+2);
        He(i+1) = HE;
        delstar(i+1)=H*theta(i);
    end
end

