% The size of the matrix A_j will be 2*NN. There are NN eigenvalues uniformly
% distributed in [-a,a] (a < 1)and 20 eigenvalues on the circle of radius
% absEvs(j) for j=1,...,4

rng(42);
NN = 20;
% The matrix A_j will be A_j = Q'*D*Q, where D is diagonal.
[Q,~] = qr(rand(2*NN));
I = eye(2*NN);
a = -0.9;
aux = linspace(-a,a,NN); % evs inside the unit disk, which is the target set
absEvs = [1.5, 2.5, 4 6 100]'; % absolute value of the evs outside the disk
K = length(absEvs);
% eigenvalues outside the disk
outEvs = absEvs*exp(2*pi*1i*rand(1,NN));
% setting down some variables
realEvs = zeros(K,2*NN);
A = zeros(2*NN,2*NN, K);
% There are 2*K matrix-valued functions to be tested. They have the form G_j(z)
% = A_j - zI and F_j(z) = G_(z)/prod(z-poles) for j=1,..,K. Notice that G_j
% are holomorphic, while F are meromorphic with poles "poles". They are
% saved in the homonym cells.
G = cell(K,1);
% The eigenvalues, backward errors and info returned by the contour
% integrals for each of the 2*K functions are saved in these variables.
% We use the suffix "-l" for the Loewner method 
EvsG = cell(K,1); EvsGl = cell(K,1);
residsG = G; residsGl = cell(K,1);
infoG = cell(K,1);   infoGl = cell(K,1);


%% Optional parameters for contour solvers
gam = 0; rad = 1;
opts.verbose = 0;
opts.thres = 1e-9;
opts.onlyIn = 1;
opts.ref2 = 0;
opts.GK = 0;
opts.maxRefine = 0;
opts.method = 'std';
opts.nc = 32;
opts.V = rand(2*NN,NN);
opts.M = 1;
% Set the Loewner data. We set the points sigma to be a dist (twice) the
% radius from the target set.
nLoewPoints = NN;
dist = 2;
lData.dir = rand(2*NN, 2*nLoewPoints);
theta = 2*pi*(1:2*nLoewPoints)/(2*nLoewPoints);
lData.points = gam + dist*exp(1i*theta);
opts.lData = lData;
%% Cycle of contour solvers
for j = 1:K
   % The real eigenvalues
   realEvs(j,:) = [aux outEvs(j,:)];
   % The matrices A_j
   A(:,:,j) = Q*diag(realEvs(j,:))*Q';
   G{j} = @(z) A(:,:,j) -z*I;
  
   opts.method = 'std';
   [EvsG{j}, EvecsG{j}, residsG{j}, infoG{j}] = contourSolver(G{j}, gam, rad, opts);
   opts.method = 'loewner';
   [EvsGl{j}, ~, residsGl{j}, infoGl{j}] = contourSolver(G{j}, gam, rad, opts);
   fprintf('\n***** Absolute values of outside evs: %2.2f *****\n', absEvs(j))
   fprintf('Std       G: # of evs: %d, max(resids): %2.2e\n', length(EvsG{j}), max(residsG{j}))
   fprintf('Loewner   G: # of evs: %d, max(resids): %2.2e\n', length(EvsGl{j}), max(residsGl{j}))   
end

%% Options for the filter functions
% Build the filter functions: b0 is the standard for hankel, bs is the
% standard for Loewner, while br and and brs are the corresponding rational
% ones
sigma = dist;
N = opts.nc;
bs = @(x) 1./(sigma-x).*((1-x.^N).^-1-(1-sigma^N).^-1);
p=0;
b0 = @(x) x.^p.*(1-(x).^N).^-1;
a = -4;
X = linspace(-a,a,1007);
%% Plotting
% All the following lines take care of drawing plots
subplot(2,2,1)
hold off
colors = lines(6);
grid
for j=1:K
    semilogy(real(EvsG{j}),residsG{j}, 'o', 'Color', colors(j,:))
    hold on
    Legend{j} = join(["$|r_", string(j), "| = ", string(absEvs(j)), " $" ]);
end
title({'$\eta(\lambda_k, v_k)$ for the Hankel algorithm'}, 'Interpreter', 'latex')
legend(Legend, 'Interpreter', 'latex', 'NumColumns', 2, 'FontSize',11, 'Location', 'best')
xlabel({'Eigenvalues $\lambda_k$'},'Interpreter', 'latex')
grid
axis([-1,1, 1e-17,1.2])
subplot(2,2,3)
hold off
grid
for j = 1:K
    semilogy(real(EvsGl{j}),residsGl{j}, 'o', 'Color', colors(j,:))
    hold on
    Legendl{j} = join(["$|r_", string(j), "| = ", string(absEvs(j)), " $" ]);
end
title({'$\eta(\lambda_k, v_k)$ for the Loewner algorithm'}, 'Interpreter', 'latex')
legend(Legendl, 'Interpreter', 'latex', 'NumColumns', 2, 'FontSize',11, 'Location', 'best')
xlabel({'Eigenvalues $\lambda_k$'},'Interpreter', 'latex')
axis([-1,1, 1e-17,1.2])
%axis([-1,1, 1e-20,1.2])

grid

subplot(2,2,[2 4])
hold off
semilogy(X,abs(b0(X)), '-', 'Color', colors(1,:))
hold on
semilogy(X,abs(bs(X)), '-', 'Color', colors(2,:))
grid
legend({'$|b_0(z)|$', '$|b_\sigma(z)|$'}, 'Interpreter', 'latex', 'Location', 'best')
axis([-4,4, 1e-20,35])


mypdf('filterMero', 0.65,2.5)

