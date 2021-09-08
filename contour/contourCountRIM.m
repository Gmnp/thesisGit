function nevs = contourCountRIM(F, n, nc, weights, lx, z, lupoints, Fp, Fpexists, verbose) 

    
    % nevs  = CONTOURCOUNT3(F, n, nc, gam, rad, w, z, lupoints, Fp,
    % ... Fpexists, verbose) is an auxiliary function of
    % CONTOURSOLVER2. It  computes an approximation of the number
    % of eigenvalues minus the number of poles of a meromorphic
    % matrix function F(z) inside a circular region Gamma. if F(z)
    % is holomorphic, then CONTOURCOUNT3 returns an approximation
    % of the number of eigenvalues. 
    %
    % The script is based on the formula
    %
    %   nevs = \frac{1}{2pi*i}*\int_Gamma trace(F(z)^{-1}F(z)') dz
    %
    % The trace operator is stochastically estimated, as described in
    % [CITE]. 
    %
    % INPUT
    % - F: the square matrix function. It is a function_handle.
    % - n: size of F
    % - nc: number of quadrature points.
    % - gam: center of Gamma.
    % - rad: radius of Gamma.
    % - w: vector of length nc containining the quadrature points.
    % - z: vector of length nc containinig the points where the
    % function F is evaluated. 
    % - lupoints: a structure of length nc with the LU
    % factorization of F(z_j) for j= 1, ..., nc. See LUFACT.
    % - Fp: If Fpexists = 1, then Fp is the function_handle with the
    % derivative of F. Otherwise, an approximation is computed during
    % runtime 
    % - Fpexists: A flag that tells whether to use Fp.
    % - verbose: a flag to increase the verbosity. It ranges from 0
    % to 2.
    %
    % OUTPUT - nevs: the number of eigenvalue minus the poles of the
    %function F(z).

    
  
    % number of random vectors
    L = max(30, n/10);
    % scaling the direction of the derivative
    delta = 1e-4;
    % If the matrix is small, then we compute the trace exactly by the
    % formula trace(A) = sum_{i=1}^n(e_i^T*A*e_i)
    if n <= 100
        Stovecs = eye(n);
        exact = 1;
    else
        Stovecs = sign(rand(n,L) - .5);
        exact = 0;
    end
    nr = 0;
    if ~Fpexists
        for k = 1:nc
            v = delta*exp(2*pi*1i*rand(1)); % We take a random
            % direction
            Fpzk = (F(z(k) + v) - F(z(k)))/v; %secant approximation
            estTrace = 0;
            for j = 1:size(Stovecs,2)
                aux = Stovecs(:,j)'/lupoints{k};
                estTrace = estTrace + aux*(Fpzk*Stovecs(:,j));
            end
            nr = nr + weights(k)*estTrace;
        end
    else
        for k = 1:nc
            estTrace = 0;
            for j = 1:size(Stovecs,2)
                aux =  Stovecs(:,j)'/lupoints{k};
                estTrace = estTrace + aux * (Fp(z(k))*Stovecs(:,j));
            end
            nr = nr + weights(k)*estTrace;
        end
    end
    if ~exact
        nr = nr/L;
    end
    nevs = 4*nr*lx/(2i*pi*nc); % we divide by the weight. Notice that nc here is 4 times nc in contourRIM
    if (verbose >= 1)
        if exact
         fprintf(['Estimation of #(eigenvalues - '...
            'poles) with %d integration points: %d\n'], ...
            nc, round(real(nevs)));    
        else
        fprintf(['Stochastic estimation of #(eigenvalues - '...
            'poles) with %d integration points and ' ...
            '%d vectors: %d\n'], ...
            nc, L, round(real(nevs)));
        end
    end
    %% number of (eigenvalues - poles)
end


