%Week 2 Exercise 2

close all
clear

ReL = [1e06, 10e06, 100e06]; %% three dif 
x = linspace(0,1,101);

uegrad = [-0.2, 0, 0.2]; %% three dif 
ue = zeros(1,length(x),length(uegrad));
ue(1,:) = 1;

for i = 2:length(x)
    for j = 1:length(uegrad)
        ue(i,j) = ue(i-1,j) + uegrad(j)*x(2);
    end
end

figure 
plot(x,ue(:,1),x,ue(:,2),x,ue(:,3))

tableX = zeros(3,3);
tableRethet = zeros(3,3);

for r = 1:length(ReL)
    for j = 1:length(uegrad)
        
        integral = 0; 
        Mthick = 0;

        n = 101;
        laminar = true;
        i = 1;
        
        while laminar && i < n
            i = i+1;
            integral = integral + ueintbit(x(i-1),ue(i-1,j),x(i),ue(i,j)); 
            Mthick = sqrt(0.45*(ue(i,j)^-6)*integral/ReL(r));
            Rethet = ReL(r)*ue(i,j)*Mthick;
           
            m = -ReL(r)*(Mthick^2)*uegrad(j);

            H = thwaites_lookup(m);
            He = laminar_He(H);

            if log(Rethet) >= 18.4*He - 21.74
                laminar = false;
                disp([x(i) Rethet/1000])
                tableX(r,j) = x(i)
                tableRethet(r,j) = Rethet/1000
            end
 
        end
        if i == n 
            disp('No Transition')
            tableX(r,j) = NaN
            tableRethet(r,j)= NaN 
        end
    end 
end 

table(:,1) = tableX(:,1); table(:,3) = tableX(:,2); table(:,5) = tableX(:,3);
table(:,2) = tableRethet(:,1); table(:,4) = tableRethet(:,2); table(:,6) = tableRethet(:,3);

