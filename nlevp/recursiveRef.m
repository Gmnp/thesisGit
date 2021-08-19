% recursive refinement
rng(0);
n_clusters = 3;
centres= linspace(-0.8, 0.8, n_clusters)*1i;
radius = 0.2;
n_b_evs = 2; % bad eigenvalues per cluster
b_evs = zeros(n_clusters, n_b_evs);
for j=1:n_clusters
    b_evs(j,:) = exp(2*pi*1i*rand(n_b_evs,1)).*rand(n_b_evs,1)*radius +centres(j);
end

n_g_evs = 4;
g_evs = exp(2*pi*1i*rand(n_g_evs,1))*1.5;

% figure(1)
% hold off
% plot(g_evs, '*')
% hold on
% plot(b_evs(:), '.')

% Build the matrix

P = cell(n_clusters+2,1);
for j=1:n_clusters
   P{j} =  b_evs(j,:);
end
P{n_clusters+1} = g_evs;
P{n_clusters+2} = centres;
[Q,~] = qr(rand(n_clusters+2));

F = @(z) Q*diag( [cellfun(@(y) prod(z-y), P(1:end-1)); prod(z-P{end})^-2 ])*Q';

clear opts;
opts.verbose = 2;
opts.V = rand(size(F(1)))-0.5;
opts.V = opts.V./vecnorm(opts.V);
opts.thres = 1e-15;
opts.ref2 = 0;
opts.GK = 0;
opts.maxRefine = 0;
opts.M = 4;
opts.nc = 16;

gam = 0; rad = 2;
[Evs, Evecs, resids, info] = contourSolver(F, gam, rad, opts);
opts.ref2 = 3;
[Evsr, Evecs, residsr, infor] = contourSolver(F, gam, rad, opts);

% Compute the forward error of the eigenvalues



%% Plots

% Eigenvalues
Evs = sort(Evs);
Evsr = sort(Evsr);
figure(1)
hold off
plot(g_evs, '*')
hold on
plot(b_evs(:), '*')
nc = 100;
w = exp(1i*2*pi*(0:nc)/nc);
plot(w*2, '-k')
legend({'Good eigenvalues', 'Bad eigenvalues', '$\partial\Omega$'}, 'Interpreter', 'latex')
mypdf('revursiveEvs', 0.7,2)


% Forward error
realEvs = sort([g_evs(:); b_evs(:)]);
f_err = abs(realEvs(1:n_clusters*n_b_evs) - Evs(1:n_clusters*n_b_evs)).*abs(realEvs(1:n_clusters*n_b_evs));
f_err_rec =  abs(realEvs(1:n_clusters*n_b_evs)- Evsr(1:n_clusters*n_b_evs)).*abs(realEvs(1:n_clusters*n_b_evs));
figure(2)
hold off
semilogy(f_err, '*');
hold on
semilogy(f_err_rec, 'o')
legend({'No refiment', 'Recursive refinement'}, 'Interpreter', 'latex', 'Location', 'southeast')
title({'Relative Forward Error'}, 'Interpreter', 'latex')
xlabel({'Index of the bad eigenvalues'}, 'Interpreter', 'latex')
mypdf('revursiveError', 0.7,2)


% Elbow plot
figure(3)
n_clusters = 5;
n_b_evs = 10; % bad eigenvalues per cluster
b_evs = zeros(n_b_evs, n_clusters);
centres= linspace(-0.8, 0.8, n_clusters);
centres = (rand(n_clusters) +1i*rand(n_clusters))*3;
for j=1:n_clusters
    b_evs(:,j) = exp(2*pi*1i*rand(n_b_evs,1)).*rand(n_b_evs,1)*0.1 +centres(j);
end
points = [real(b_evs(:)), imag(b_evs(:))];
for k = 2:8
   [~,~,er] = kmeans(points,k);
   ers(k) = sum(er);
end

semilogy([2:8], ers(2:8), '-o')
title({'The elbow plot'}, 'Interpreter', 'latex')
xlabel({'Number of clusters k'}, 'Interpreter', 'latex')
legend({'$\Theta(k)$'}, 'Interpreter', 'latex')
mypdf('elbow', 0.7,2)


figure(4)
hold off
plot(b_evs, 'o')
axis([0.5, 2.5, -0.2, 2.7])
mypdf('clusters', 0.7,2)
