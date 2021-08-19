% This script prints the tropical function with built with the harmonics
% and its Newton polygon.

nTrunc = 1000; % where we truncate the \tropx f(x)
nRoots = 25; % number of tropical roots and number of Newton vertices
harmon= cumsum([1:nRoots].^-1)';
coeffs = exp(cumsum([1:nTrunc].^-1)');

% Newton polygon
figure(1)
hold off
plot(harmon, '-*')
%grid
%legend({'Newton polygon'}, 'Interpreter', 'latex', 'Location', 'NorthWest')

nPoints = 10000;
x = linspace(0, 0.999, nPoints);
% This is the definition of the tropical function
f = @(z) max(coeffs.*(z.^([1:nTrunc]'*ones(1,length(z)))));
fx = f(x);
%the tropical roots
troots = [0 exp(-[1:nRoots].^-1)];

% Values of the function \tropx f(x)
figure(2)
hold off
plot(x,fx, '-')
gtroots = f(troots);
hold on
plot(troots, gtroots, '.', 'MarkerSize',30);
%grid
legend({'$\mathsf{t}_{\times} f(x)$', 'Tropical roots'}, 'Interpreter', 'latex', 'Location', 'NorthWest')
axis([0, 1, 0, 16])
% Zoom of values of the function \tropx f(x)
figure(3)
hold off
plot(x,fx, '-')
gtroots = f(troots);
hold on
plot(troots, gtroots, '.', 'MarkerSize',30);
%grid
legend({'$\mathsf{t}_{\times} f(x)$', 'Tropical roots'}, 'Interpreter', 'latex', 'Location', 'NorthWest')
axis([0.6, 0.87, 1.6, 5])