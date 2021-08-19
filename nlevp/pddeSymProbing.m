clear opts
opts.verbose=2;
opts.ref2 = 0;
opts.nc = 32;
opts.maxRefine = 0;
opts.M = 1;
rng(42)
P1 = rand(81,25);
P1 = P1./vecnorm(P1);
P2 =  rand(81,12);
P2 = P2./vecnorm(P2);
opts.V = P1;

% Under these settings, photonic has 28 eigenvalues.
gam = 0; rad = 1.5;
[coeffs,fun,F] = nlevp( 'pdde_symmetric',10);

[Evs, EvecsA, residsrng1, info] = contourSolver(F, gam,rad,opts);
opts.M = 2;
[Evs, Evecs, residsrng2, info] = contourSolver(F, gam,rad,opts);
opts.M = 1;
opts.V = EvecsA;
[Evs, Evecs, residsEvec1, info] = contourSolver(F, gam,rad,opts);
opts.M = 2;
[Evs, Evecs, residsEvec2, info] = contourSolver(F, gam,rad,opts);
opts.M = 1;
opts.V = [P1, P2];
[Evs, Evecs, residsrng3, info] = contourSolver(F, gam,rad,opts);
opts.V = [EvecsA, P2];
[Evs, Evecs, residsEvec3, info] = contourSolver(F, gam,rad,opts);

%% Plots
colors = lines(6);
figure(5)
hold off
semilogy(residsrng1, '--x', 'Color', colors(1,:))
hold on
semilogy(residsrng2, '--*', 'Color', colors(1,:) )
semilogy(residsEvec1, '--x', 'Color', colors(2,:))
semilogy(residsEvec2, '--*', 'Color', colors(2,:))
semilogy(residsrng3, '--x', 'Color', colors(6,:))
semilogy(residsEvec3, '--x', 'Color', colors(4,:))
grid
legend({'$P = P_1$, $m=1$', '$P = P_1$, $m=2$', '$P = V$, $m=1$', '$P = V$, $m=2$', '$P = [P_1\ P_2]$, $m=1$', '$P = [V\ P_2]$, $m=1$'},'NumColumns', 3, 'Interpreter', 'latex', 'Location', 'south')
axis([0,25, 1e-12,1e-1])
mypdf('pdde32', 0.7,1.8)


%% nc= 128
opts.nc = 64;

opts.V = P1;
[Evs, EvecsA, residsrng1, info] = contourSolver(F, gam,rad,opts);
opts.M = 2;
[Evs, Evecs, residsrng2, info] = contourSolver(F, gam,rad,opts);
opts.M = 1;
opts.V = EvecsA;
[Evs, Evecs, residsEvec1, info] = contourSolver(F, gam,rad,opts);
opts.M = 2;
[Evs, Evecs, residsEvec2, info] = contourSolver(F, gam,rad,opts);
opts.M = 1;
opts.V = [P1, P2];
[Evs, Evecs, residsrng3, info] = contourSolver(F, gam,rad,opts);
opts.V = [EvecsA, P2];
[Evs, Evecs, residsEvec3, info] = contourSolver(F, gam,rad,opts);


colors = lines(6);
figure(6)
hold off
semilogy(residsrng1, '--x', 'Color', colors(1,:))
hold on
semilogy(residsrng2, '--*', 'Color', colors(1,:) )
semilogy(residsEvec1, '--x', 'Color', colors(2,:))
semilogy(residsEvec2, '--*', 'Color', colors(2,:))
semilogy(residsrng3, '--x', 'Color', colors(6,:))
semilogy(residsEvec3, '--x', 'Color', colors(4,:))
grid
legend({'$P = P_1$, $m=1$', '$P = P_1$, $m=2$', '$P = V$, $m=1$', '$P = V$, $m=2$', '$P = [P_1\ P_2]$, $m=1$', '$P = [V\ P_2]$, $m=1$'},'NumColumns', 3, 'Interpreter', 'latex', 'Location', 'best')
axis([0,25, 1e-12,1e-1])
mypdf('pdde128', 0.7,1.8)

