% Script to produce a very easy toy example
n = -4; m = 4;
jn = [n:-1];
jp = [0:m];
bjn = 2*3.^(jn+2);
bjp = 3*2.^(-jp);

plot([jn jp], log([bjn bjp]), '-*')
%grid
%legend({'Newton Polygon'},'Interpreter','latex')
