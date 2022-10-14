%% Solar collector
a=1;
b=1;
c=1;

z=-1:0.01:1;
x=a*(z-b).^2+c;


plot(x,z)
xlabel('x')
ylabel('z')
title('$x=a(z-b)^2+c$', 'interpreter','latex')
% syms x y
% fimplicit(0.211*x^2-0.157*x*y+0.065*y^2==1, [-4 4 -6 6])