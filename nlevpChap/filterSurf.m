a = -3;
nn = 507;
X = linspace(a,-a,nn);
Y = linspace(-a*1i, a*1i, nn);
Z = X+Y.';

%% Define the filter function
N = 32;
sigma  = -2;
bs = @(x) 1./(sigma-x).*((1-x.^N).^-1-(1-sigma^N).^-1);
b0 = @(x) (1-(x).^N).^-1;

%% Surfaces scaled logarithmically
Z0 = (abs(b0(Z)));
Zs = (abs(bs(Z)));

%% Plots
colors = lines(6);

figure(1)
hold off
surf(X,imag(Y.'), Z0, 'FaceColor',colors(1,:), 'FaceAlpha',0.7, 'Edgecolor', 'None')
hold on
surf(X,imag(Y.'), Zs, 'FaceColor',colors(2,:), 'FaceAlpha',0.7, 'Edgecolor', 'None')
legend({'$|b_0(z)|$', '$|b_\sigma(z)|$'}, 'Interpreter', 'latex')
set(gca,'ZScale','log')
mypdf('surfFilter1', 0.6,1.7)

figure(2)
hold off
semilogy(X,abs(b0(X)), '-', 'Color', colors(1,:))
hold on
plot(X,abs(bs(X)), '-', 'Color', colors(2,:))
legend({'$|b_0(z)|$', '$|b_\sigma(z)|$'}, 'Interpreter', 'latex')
grid
mypdf('surfFilter2', 0.6, 1.7)
