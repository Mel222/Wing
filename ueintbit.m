function f = ueintbit(xa,ua,xb,ub)

u_bar = (ua+ub)/2;
delta_u = ub - ua;
delta_x = xb - xa;

f = (u_bar^5+(5/6)*(u_bar)^3*(delta_u)^2+(1/16)*u_bar*(delta_u)^4)*delta_x;