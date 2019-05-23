clear 

close all



global Re ue0 duedx

ReRange= [1e6, 1e7, 1e8,1e7,1e7];

ue0=1;  %free stream velocity

duedxRange = [-0.6, -0.6, -0.6,-0.3,-0.9];

Comb = [duedxRange;ReRange];

TableSep = zeros(1,length(Comb(1,:)));



figure(1)



for c = 1:length(Comb(1,:))

    Re = Comb(2,c)

    duedx = Comb(1,c)



    %Definition of the initial value of theta and deltaE

    x0=0.01;

    thick0(1)=0.037*x0*(Re*x0)^(-1/5);

    thick0(2)=1.80*thick0(1);



    %the differential equation is solved



    [delx, thickhist]=ode45(@thickdash,[0 0.99], thick0);



    x=x0+delx;



    HE=thickhist(:,2)./thickhist(:,1);



    seplocation=0;

    i=1;



    while i<length(HE) && seplocation==0

        if HE(i)<=1.46

            seplocation=x(i);

        end

        i=i+1;

    end  



    hold on 

    plot(x,HE);



    TableSep(c) = seplocation;

    

    if Re == 1e7 && duedx == -0.6

        Thetacurve = thickhist(:,1)

        Ethickcurve = thickhist(:,2)
        
        xcurve = x

        

    end

    

end



legend(num2str(Comb(:,1)),num2str(Comb(:,2)),num2str(Comb(:,3)),num2str(Comb(:,4)),num2str(Comb(:,5)))

hold off 



figure(2)

plot(xcurve, Thetacurve,xcurve,Ethickcurve)

legend('theta','deltaE')

