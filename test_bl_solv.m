%  Zero pressure gradient b-l; bl_solv and Blasius

clear all
global Re

% Re = 10000000;
% Re = 1e6;
Re = 2500;

n = 100;
x = linspace(1/n,1,n);
cp = zeros(1,n);

[int, ils, itr, its, delstar, theta] = bl_solv ( x, cp );
disp([num2str(int), num2str(ils),num2str(itr),num2str(its)])
blthet = 0.664 * sqrt(x/Re);

if int~=0
  disp(['Natural transition at x = ' num2str(x(int))])
end

plot(x,blthet,'-',x,theta,'x')
xlabel('x')
ylabel('\theta')
legend('Blasius','blsolv')