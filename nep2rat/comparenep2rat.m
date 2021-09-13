% We want to compare some algorithms:
% - weighted  AAA
% - surrogate AAA + Leja-Bagby
% - contourHankel
% - contourLoewner
% - contour* + refinements (does not matter HankelLoewner)
% - maybe GK method?

% We will use the same problems of the previous paper

% clamped and photonic are slow

clear optsnep2rat optsCS

%% Extract all reps and neps from NLEVP.
pep = nlevp('query','pep');
allProblems = nlevp('query','problems');
nep = setdiff(allProblems, pep);
nep = setdiff(nep, {'fiber', 'gun', 'pillbox_cavity', 'sandwich_beam', 'laser'}); 
%nep = {'schrodinger_abc'}
nb_test_pbs = length(nep);
fprintf('Number of test problems: %3.0f\n',nb_test_pbs)

%% General settings
N = 10;
nRows = ceil(sqrt(nb_test_pbs));
nCols = ceil(nb_test_pbs/nRows);
nameAlgs = {'Weighted', 'AAA+LB', 'Hankel', 'Loewner', 'Refined'};
onlyIn = 1;

%% General settings for nep2rat
tol = 1e-7;
nc = 300; %number of sample points for Sigma
nc2 = 50; %number of sample points on the countour of Sigma
dmax = 100; %max degree
evsThresh = 1e-3; % Removing the bad eigenvalues
% Option parameters for nep2rat.m
useZ2 = 1;
optsnep2rat.dmax = dmax;
optsnep2rat.tol1 = tol;
optsnep2rat.tol2 = tol;
optsnep2rat.evsThresh = evsThresh;
optsnep2rat.verbose = 0;

%% General settings for contourSolver
optsCS.GK = 0;
optsCS.onlyIN = onlyIn;
optsCS.verbose = 0;

%% Matrices of output and other memory allocations
nAlg = 5;
timings = zeros(nb_test_pbs, nAlg);
residss = cell(nb_test_pbs, nAlg);
Evss = residss;
Evecss = residss;
infos = residss;
gam = zeros(nb_test_pbs);
rad = gam;
%% Cycle
for kk = 1:nb_test_pbs
    
    switch nep{kk} %generate F
        
        case 'bent_beam' % temporary
            gam(kk) = 60;
            rad(kk) = 30;
            %half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';
            
        case 'buckling_plate'  % The smallest poles are in k*pi/2, and in 4.70
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';
        
             case 'canyon_particle'
            gam(kk) = -9e-2+1e-6i;
            rad(kk) = .1;
            %half_disc(kk) = 1; %half disc domain
            stepSize = 1;
            [coeffs,fun,F] = nlevp(nep{kk}, stepSize);
             optsnep2rat.nevs = 'all';
    
            
        case 'clamped_beam_1d'
            gam(kk) = 0;
            rad(kk) = 10;
            [coeffs,fun,F] = nlevp(nep{kk}, 100);
            optsnep2rat.nevs = 'all';

      case 'distributed_delay1'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'fiber'
            gam(kk) = 0;
            rad(kk) = 2e-3;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 2;

      case 'gun'
            gam(kk) = 62500;
            rad(kk) = 50000;
            half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 30;
          

      case 'hadeler'
            gam(kk) = -30;
            rad(kk) = 11.5;
            [coeffs,fun,F] = nlevp(nep{kk},200);
            optsnep2rat.nevs = 'all';

      case 'loaded_string'
            gam(kk) = 362;%14;
            rad(kk) = 358;%11;
            %KAPPA = 1; mass = 1; % pole is KAPPA/mass
            [coeffs,fun,F] = nlevp(nep{kk},100);%, N, KAPPA, mass);
            optsnep2rat.nevs = 'all';

      case 'nep1'
            gam(kk) = 0;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'nep2'
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case  'nep3'
            gam(kk) = 5i;
            rad(kk) = 2; % found 14 evals in this disc
            [coeffs,fun,F] = nlevp(nep{kk},N);
            optsnep2rat.nevs = 'all';

      case 'neuron_dde'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'photonic_crystal'
            % The poles of the default values are +-1.18-0.005i and +-1.26-0.01
            gam(kk) = 11;
            rad(kk) = 9;
            [coeffs,fun,F] = nlevp(nep{kk}, N);
            optsnep2rat.nevs = 'all';

            
      case 'pillbox_small'
            gam(kk) = 0.08;
            rad(kk) = 0.05;
            %half_disc(kk) = 1; %half disc domain
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'sandwich_beam'
            gam(kk) = 7000;
            rad(kk) = 6900;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 200;
       
      case 'railtrack2_rep'
            % railtrack2_rep has a pole in 0, so we shifted it in 5 (?)
            gam(kk) = 3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk}, N);
            optsnep2rat.nevs = 'all';

      case 'railtrack_rep'
            gam(kk) = -3;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'schrodinger_abc'
            gam(kk) = -10;
            rad(kk) = 5;
            [coeffs,fun,F] = nlevp(nep{kk}, N);
             optsnep2rat.nevs = 'all';
           
      case 'square_root'
            gam(kk) = 10+50i;
            rad(kk) = 50;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'time_delay'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'time_delay2'
            gam(kk) = 0;
            rad(kk) = 15;
            [coeffs,fun,F] = nlevp(nep{kk});
            optsnep2rat.nevs = 'all';

      case 'time_delay3'
            gam(kk) = 2;
            rad(kk) = 3;
            [coeffs,fun,F] = nlevp(nep{kk}, N, 5);
            optsnep2rat.nevs = 'all';

      otherwise
            gam(kk) = 0;
            rad(kk) = 2;
            [coeffs,fun,F] = nlevp(nep{kk}, N);
            optsnep2rat.nevs = 'all';

    end
    
    Fsize(kk) = length(coeffs{1}); %record the size of each NEP
    %fprintf('Pb size %5d\n', Fsize(kk))
    fprintf('*******************************\n');
    fprintf('kk=%2d, Problem: %s, n = %4d\n',kk, nep{kk},Fsize(kk));
    fprintf('*******************************\n');
    
    % Generate ZZ and Z2
    rng(0); %Fix the random number generator
    Z(kk,:) = rand(1,nc).*exp(rand(1,nc)*2*pi*1i);
    Z2 = gam(kk) + rad(kk)*exp(1i*linspace(0,2*pi,2*nc2));
    Z(kk,:) = Z(kk,:)*rad(kk)  + gam(kk); % shift to the correct points
    
     if useZ2
       optsnep2rat.Z2 = Z2;
     end
    
     alg = 0;
    
     %% Weighted
     disp('Weighted')
     alg = alg+1;
     Fw.fun = fun;
     Fw.coeffs = coeffs;
     optsnep2rat.phase1 = 'weighted';
     optsnep2rat.phase2 = '';
     
     tstart = tic;
     [Evs, Evecs, resids, info] = mixedSolver(Fw, Z(kk,:), optsnep2rat);
     timings(kk,alg) = toc(tstart);
     if onlyIn
         flag = isInsideContour(Evs, gam(kk), rad(kk));
         Evs = Evs(flag);
         resids = resids(flag);
         Evecs = Evecs(:,flag);
     end
     residss{kk,alg} = computeResids(F, Fw, Evs, Evecs);
     Evss{kk,alg} = Evs;
     Evecss{kk,alg} = Evecs;
     infos{kk,alg} = info;
     
     
     %% Surrogate plus Leja-Bagby
     disp('Surrogate AAA+Exact')
     alg = alg+1;
     optsnep2rat.phase1 = 'sur';
     optsnep2rat.phase2 = 'LB';
     
     tstart = tic;
     [Evs, Evecs, resids, info] = mixedSolver(F, Z(kk,:), optsnep2rat);
     timings(kk,alg) = toc(tstart);
     if onlyIn
         flag = isInsideContour(Evs, gam(kk), rad(kk));
         Evs = Evs(flag);
         resids = resids(flag);
         Evecs = Evecs(:,flag);
     end
     residss{kk,alg} = computeResids(F, Fw, Evs, Evecs);
     Evss{kk,alg} = Evs;
     Evecss{kk,alg} = Evecs;
     infos{kk,alg} = info;
     
     %% ContourHankel
     disp('Contour Hankel')
     alg = alg+1;
     optsCS.method = 'std';
     optsCS.ref2 = 0;
     optsCS.maxrefine = 0;
     
     tstart = tic;
     [Evs, Evecs, resids, info] = contourSolver(F, gam(kk), rad(kk), optsCS);
     timings(kk,alg) = toc(tstart);
     residss{kk,alg} = computeResids(F, Fw, Evs, Evecs);
     Evss{kk,alg} = Evs;
     Evecss{kk,alg} = Evecs;
     infos{kk,alg} = info;     
     
     
     %% ContourLoewner
     disp('Contour Loewner')
     alg = alg+1;
     optsCS.method = 'loewner';
     
     tstart = tic;
     [Evs, Evecs, resids, info] = contourSolver(F, gam(kk), rad(kk), optsCS);
     timings(kk,alg) = toc(tstart);
     residss{kk,alg} = computeResids(F, Fw, Evs, Evecs);
     Evss{kk,alg} = Evs;
     Evecss{kk,alg} = Evecs;
     infos{kk,alg} = info;
     
     
     %% ContourRefined
     disp('Contour Refined')
     alg = alg+1;
     optsCS.method = 'std';
     optsCS.ref2 = 2; % later, it will be 1
     optsCS.maxrefine = 1;
     
     tstart = tic;
     [Evs, Evecs, resids, info] = contourSolver(F, gam(kk), rad(kk), optsCS);
     timings(kk,alg) = toc(tstart);
     residss{kk,alg} = computeResids(F, Fw, Evs, Evecs);
     Evss{kk,alg} = Evs;
     Evecss{kk,alg} = Evecs;
     infos{kk,alg} = info;

end


%% Plotting
symbPlot = '+xsdo';
figure(1)
hold off
for kk = 1 : nb_test_pbs
    subplot(nRows, nCols, kk)
    for alg = 1:nAlg
        if isempty(residss{kk, alg})
            plot(NaN, symbPlot(alg));
        else
            semilogy(residss{kk,alg}, symbPlot(alg))
        end
        hold on
    end
    grid
%    legend([nameAlgs], 'Interpreter', 'latex')
    myTitle = replace(nep{kk}, '_', '\_');
    title([myTitle], 'Interpreter', 'latex')
end
    legend([nameAlgs], 'Interpreter', 'latex', 'Location', 'best')

    % timings
figure(2)
hold off
for alg = 1:nAlg
        semilogy(timings(:, alg), strcat('--', symbPlot(alg)))
        hold on
 end
grid
legend([nameAlgs], 'Interpreter', 'latex', 'Location', 'best')
title('Timings of the problems', 'Interpreter', 'latex')
axis([0, 22, 10^-2,10^4])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
nepR =  cellfun(@(x)replace(x, '_', '\_'), nep, 'UniformOutput', false);
xticks([1:21])
xticklabels(nepR) 
xtickangle(90)

%% Auxiliary functions
function normFSigma = computeNorm(F,Z)
% computeNorm computes an approximation normFSigma of ||F||_Z := max_{z\in Z}
% ||F(z)||_2 and return the matrix-valued. We assume F is a struct.
    sparseFlag = issparse(F.coeffs{1});
    n = size(F.coeffs{1},1);
    d = length(F.coeffs);
    state = rng();
    rng('default');
    u = randn(n,1);
    u = u/norm(u);
    rng(state);
    if sparseFlag
        uj = sparse(n,d);
    else
        uj = zeros(n,d);
    end
    for jj = 1:d
        uj(:,jj) = F.coeffs{jj}*u;
    end
    normFSigma = 0;
    for jj =1:length(Z)
        normFSigma = max(normFSigma, norm(uj*F.fun(Z(jj)).'));
    end
end

function resids = computeResids(F, Fw, Evs, Evecs)
n = length(Evs);
resids = zeros(1,n);
alphas = cellfun(@(z) norm(z, 'fro'), Fw.coeffs);
for j = 1:n
  resids(j) = norm(F(Evs(j))*Evecs(:,j))/(alphas*abs(Fw.fun(Evs(j)))');
end 
end
