% Script to produce a very easy toy example: Newton polygon of the
% polynomial [1e-70 0 0 0 1e39 1e3 1e10 1e40] plus the holomorphic function
% in A(1,a) (1-a)/((z-1)(z-a)), where a = 1e28.

n = -4; m = 4;
jn = [n:-1];
jp = [0:m];
bjn = 2*3.^(jn+2);
bjp = 3*2.^(-jp);

plot([jn jp], log([bjn bjp]), '-*')
%grid
%legend({'Newton Polygon'},'Interpreter','latex')
