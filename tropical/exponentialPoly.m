% This code compute produces an example of a quartic polynomial plus an
% exponential

% construction of the polynomial
rng(42); n = 20;
B0 = randn(n); B1 = 1e5*randn(n);
B2 = 1e1*randn(n); B3 = 1e-2*randn(n); B4 = 1e3*randn(n);
% and of the exponential
C = rand(n);
I = eye(n);
P = @(z) B0 +z*(B1 + z*(B2 + z*(B3+ z*B4)));
F = @(z) z*I +exp(-z)*C;
H = @(z) F(z) + P(z);


% Express H(z) as laurent series
% Length of Truncation
k = 22;
% C1 is a vector containing the norm of the coefficient matrices of H(z)
C1 = norm(C)./factorial([0:k]);
C1(1) = C1(1) + norm(B0);
C1(2) = C1(2) + norm(B1) + 1;
C1(3) = C1(3) + norm(B2);
C1(4) = C1(4) + norm(B3);
C1(5) = C1(5) + norm(B4);

% Plot the Newton polygon
figure(1)
colors = lines(6);
hold off
plot(0:k, log(C1), '*', 'Color', colors(1,:))
aux = sort(convhull([0:k], log(C1)));
hold on
plot(aux(2:end)-1, log(C1(aux(2:end))), '*-', 'Color', colors(2,:))
% plot(0:4, log([norm(B0), norm(B1), norm(B2), norm(B3), norm(B4)]), 'd', 'Color', colors(3,:))
% plot([0 1 4], log([norm(B0), norm(B1), norm(B4)]), 'x-', 'Color', colors(4,:))
legend({'$(j, \log(||C_j||))$', 'Newton Polygon'}, 'Interpreter', 'latex', 'Location', 'best')
axis([0 k -50 15])
%grid

% Check that this example works fine

% Tropical roots for P(z)
a1P = (norm(B0)/norm(B1));
a2P = (norm(B1)/norm(B4))^(1/3);
fprintf('The first tropical roots of H(z) are %e and %f.\n', a1P,a2P)
if a1P/a2P < (1+2*cond(B0))^-2
   fprintf('The condition on the first tropical root of P(z) is satisfied.\n\n')
   radP = (1+2*cond(B0))*a1P;
   outerRadP = a2P/(1+2*cond(B0));
   fprintf('No evs in annulus A(%e,%e).\n\n', radP, outerRadP)
end

% Tropical roots for H(z)
a1 = (C1(1)/C1(2));
a2 = (C1(2)/C1(5))^(1/3);
a3 = (C1(5)/C1(22))^(1/17);
fprintf('The first tropical roots of H(z) are %e, %f, and %f.\n', a1,a2,a3)
if a1/a2 < (1+2*cond(B0 + C))^-2
   fprintf('The condition on the first tropical of H(z) root is satisfied.\n\n')
   rad = (1+2*cond(B0 + C))*a1;
   outerRad = a2/(1+2*cond(B0 + C));
   fprintf('No evs in annulus A(%e,%e).\n\n', rad, outerRad)
end






% Contour_solver is needed for this part of the code, which compute the
% eigenvalues of P and H. It can be commented out.

opts.verbose=0;
opts.nc=10;
opts.M=1;
opts.V = rand(20,20);
opts.V = opts.V./vecnorm(opts.V);
opts.ref2=0;
opts.maxRefine=0;
[Evs, Evecs, resids, info] = contourSolver(H, 0, rad, opts);

opts.nc = 18;
[EvsP, EvecsP, residsP, infoP] = contourSolver(P, 0, radP, opts);

