function [cl, cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)

cl = -2*circ;

ueTE = (1-cp(end))^0.5;
MthickTE = thetal(end)+thetau(end);
H_TE = delstarl(end)/thetal(end) + delstaru(end)/thetau(end);

MthickInf = MthickTE * (ueTE)^((H_TE+5)/2);

cd = 2*MthickInf
