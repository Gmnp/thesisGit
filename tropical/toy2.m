% Script to produce a very easy toy example: Newton polygon of the
% polynomial [1e-70 0 0 0 1e39 1e3 1e10 1e40] plus the holomorphic function
% in A(1,a) (1-a)/((z-1)(z-a)), where a = 1e28.

lx = horzcat([0 0 0 40 10 3 39], [-140:-28:-476]);
lx(11) = -70;
x = [-3:16];
hold off
plot(x,lx, '.b')
k = [4,7,11];
hold on
plot(x(k),lx(k), '*b')
plot(x(k), lx(k), '-r')
legend({'$\log_{10}$ coeffs', 'Vertices NP', 'NP'},'Interpreter','latex')
axis([-5, 17, -500,70])