
function dthickdx=thickdash(xmx0, thick)
global Re ue0 duedx
He=thick(2)/thick(1);
Rethet=Re*ue0*thick(1);
  if He>=1.46
      H=(11*He+15)/(48*He-59); 
  else
      H=2.803;
  end
cf=0.091448*((H-1)*Rethet)^(-0.232)*exp(-1.260*H);
cdiss=0.010019*((H-1)*Rethet)^(-1/6);
dthickdx=transpose([cf/2-(H+2)/ue0*duedx*thick(1), cdiss-3/ue0*duedx*thick(2)]); %thickdash must return a column vector


