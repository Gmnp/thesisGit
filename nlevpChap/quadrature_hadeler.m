clear Evs;
clear Evecs;
clear residss;
[coeffs,fun,F] = nlevp('hadeler', 30);
rng(2)
F = @(z) full(F(z));
% Setting all the parameters to find 54 evs

qPoints = [8 32  128];
nq = length(qPoints);
clear opts;
opts.verbose = 0;
opts.V = rand(size(F(1)))-0.5;
opts.V = opts.V./vecnorm(opts.V);
opts.thres = 1e-13;
opts.ref2 = 0;
opts.GK = 0;
opts.maxRefine = 0;
residss = cell(1,nq);
errQuad  = cell(1,nq);

gams = [0];
rads = [5];
Moms = [2];

nExp = 0;
%% First experiment
nExp = nExp+1;
gam = gams(nExp); rad = rads(nExp);
opts.M = Moms(nExp);
for j=1:length(qPoints)
    opts.nc = qPoints(j);
    [Evs{j}, Evecs, residss{j}, info] = contourSolver(F, gam, rad, opts);
         errQuad{j} = info.errQuad;
         
end

% We shift the color by one in the plot so we are consistent with the
% thesisButt experiment
colors = lines(6);


figure(2)
grid
% Plotting and writing legend
hold off
for j = 1:nq
semilogy(residss{j}, 'o', 'Color', colors(j,:))
hold on
Legend{j} = ['$N = ', num2str(qPoints(j)), '$.'];
end
grid
aux = axis();
axis([0 length(Evs{3})+4 1e-21 3])
%legend(Legend, 'Interpreter', 'latex')
title(['Eigenvalues backward error $\eta(\lambda_k, v_k)$'],'Interpreter', 'latex')
xlabel({'$k$'}, 'Interpreter', 'latex')
legend(Legend, 'Interpreter', 'latex', 'Location', 'best')

figure(3)
hold off
for j = 1:nq
semilogy([0:size(errQuad{j},2)-1], errQuad{j}(1,:), '-*','Color', colors(j,:))
hold on
Legend1{j} = ['$N = ', num2str(qPoints(j)), '$.'];
end
%semilogy([0:size(errQuad{nq},2)-1], errQuad{nq}(2,:), '-*')
grid
%Legend1{nq+1} = ['$N = ', num2str(2*qPoints(nq)), '$.'];
legend(Legend1, 'Interpreter', 'latex', 'Location', 'best')
title(['Relative error $E_{(N, 2N)}(A_k)$ of the quadrature'], 'Interpreter', 'latex')
xlabel({'$k$'}, 'Interpreter', 'latex')
aux = axis();
axis([0 size(errQuad{1},2)-1 aux(3:4)]);

