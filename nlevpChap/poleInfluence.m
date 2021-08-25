rng(42)
gam = 0; rad = 1;
n = 100;
xis = linspace(rad*0.5, rad*0.997, n);
% xis = 1-logspace(-1, 0,n).^-1/10000;
% xis = 1-logspace(-3, 0,n/2+1).^-1/1000;
xis = xis(2:end);
errs = zeros(4,length(xis));
errs2 = errs;
clear opts
opts.M = 2;
opts.verbose = 0;
opts.ref2 = 0;
opts.GK = 0;
opts.maxRefine = 0;
opts.V = rand(4,4)-0.5;
opts.V = opts.V./vecnorm(opts.V);
opts.nc = 6;

% We repeat the experiment N times to smooth the error
N = 1;
for k = 1:N
    for j = 1:length(xis)
        xi = xis(j);
        F4 = @(z) [z (z-xi)^-1 0 0; 0 1 (z-xi)^-1 0; 0 0 1 (z-xi)^-1; 0 0 0 1];
        F3 = @(z) [z (z-xi)^-1 0 0; 0 1 (z-xi)^-1 0; 0 0 1 0; 0 0 0 1];
        F2 = @(z) [z (z-xi)^-1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        F1 = @(z) [z 0 0 0; 0  (z-xi)^-1 0 0; 0 0 1 0; 0 0 0 1];
        [Evs4, ~, ~, info] = contourSolver(F4, gam, rad, opts);
        [Evs3, ~, ~, info] = contourSolver(F3, gam, rad, opts);
        [Evs2, ~, ~, info] = contourSolver(F2, gam, rad, opts);
        [Evs1, ~, ~, info] = contourSolver(F1, gam, rad, opts);
        errs(1,j) = errs(1,j)+min(abs(Evs4));
        errs(2,j) = errs(2,j)+min(abs(Evs3));
        errs(3,j) = errs(3,j)+min(abs(Evs2));
        errs(4,j) = errs(4,j)+min(abs(Evs1));
    end
end
errs = errs/N;

opts.nc = 64;
for k = 1:N
    k
    for j = 1:length(xis)
        xi = xis(j);
        F4 = @(z) [z (z-xi)^-1 0 0; 0 1 (z-xi)^-1 0; 0 0 1 (z-xi)^-1; 0 0 0 1];
        F3 = @(z) [z (z-xi)^-1 0 0; 0 1 (z-xi)^-1 0; 0 0 1 0; 0 0 0 1];
        F2 = @(z) [z (z-xi)^-1 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
        F1 = @(z) [z 0 0 0; 0  (z-xi)^-1 0 0; 0 0 1 0; 0 0 0 1];
        [Evs4, ~, ~, info] = contourSolver(F4, gam, rad, opts);
        [Evs3, ~, ~, info] = contourSolver(F3, gam, rad, opts);
        [Evs2, ~, ~, info] = contourSolver(F2, gam, rad, opts);
        [Evs1, ~, ~, info] = contourSolver(F1, gam, rad, opts);
        errs2(1,j) = errs2(1,j)+min(abs(Evs4));
        errs2(2,j) = errs2(2,j)+min(abs(Evs3));
        errs2(3,j) = errs2(3,j)+min(abs(Evs2));
        errs2(4,j) = errs2(4,j)+min(abs(Evs1));
    end
end
errs2 = errs2/N;

%% Plotting
figure(4)
colors = lines(6);
mSize=15;
hold off
semilogy(xis, errs(1,:), '.-', 'Color', colors(1,:), 'MarkerSize', mSize)
hold on
semilogy(xis, errs(2,:), '.-', 'Color', colors(2,:), 'MarkerSize', mSize)
semilogy(xis, errs(3,:), '.-', 'Color', colors(3,:), 'MarkerSize', mSize)
semilogy(xis, errs(4,:), '.-', 'Color', colors(4,:), 'MarkerSize', mSize)
semilogy(xis, errs2(1,:), '.--', 'Color', colors(1,:), 'MarkerSize', mSize)
semilogy(xis, errs2(2,:), '.--', 'Color', colors(2,:), 'MarkerSize', mSize)
semilogy(xis, errs2(3,:), '.--', 'Color', colors(3,:), 'MarkerSize', mSize)
semilogy(xis, errs2(4,:), '.--', 'Color', colors(4,:), 'MarkerSize', mSize)
% write legend
Legend = {'$F_4(z)$, $N=6$', '$F_3(z)$, $N=6$','$F_2(z)$, $N=6$', '$F_1(z)$, $N=6$', '$F_4(z)$, $N=64$', '$F_3(z)$, $N=64$', '$F_2(z)$, $N=64$', '$F_1(z)$, $N=64$'};
legend(Legend, 'Interpreter', 'latex', 'Location', 'best',  'NumColumns', 2)
grid
xlabel({'$\xi$'}, 'Interpreter', 'latex')
ylabel({'$|\widehat{\lambda} - \lambda|$'}, 'Interpreter', 'latex')
axis([0.5, 1, 1e-18, 0.2])
