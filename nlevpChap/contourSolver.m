function [evs, evecs, resids, info] = contourSolver(F, gam, rad, varargin)

% [Evs, Evecs, Resids, Info] = CONTOURSOLVER(F, Center, Radius) computes
% the eigenvalues and eigenvectors of a meromorphic function F(z) inside a
% disk Omega of center 'Center' and radius 'Radius'. The main idea is
% approximating a collection of matrices, called moments, defined as
% follow:
%
%        A_p = \int_\Gamma z^p*F(z)^-1*V dz,
%
% where V is a projection subspace.
%
% ============================== BASIC INPUT ==============================
%
% - F is a n-by-n matrix function_handle.
%
% - Center and Radius are the center and the radius of the contour of
% interest.
%
% =============================== OUTPUT ==================================
%
% - Evs is a vector containing the eigenvalues.
%
% - Evecs is a matrix containing  the corresponding  eigenvectors.
%
% - Resids is a vector containing an upper bound of the backward errors. We
% use the following formula: Resids(j) = norm(F(evs(j)) * evecs(:,j)) /
% norm(F)_Omega, where norm(F)_Omega = sup_{z \in Omega} norm(F(z))is
% approximated at runtime.
%
% - Info is a structure which contains information about the code. More
% details in OPTIONAL PARAMETERS.
%
% ========================== OPTIONAL PARAMETERS ==========================
%
% CONTOURSOLVER(F, Center, Radius, Opts) accepts a fourth parameter Opts.
% Opts is a structure with the following fields:
%
% .verbose [0,1,2]: Increase the verbosity of the code. Default is 0
%
% .nc [Positive even integer]: Number of trapezoidal integration points.
% Default value is 2*max(round(16*2^log10(rad)), 2).
%
% .str [positive real]: Ratio between the vertical and the horizontal axis
% of the elliptical contour. If it is given, then the horizontal axis is
% equal to "rad" and the vertical is " opts.str*rad". Default value is 1
% (circle).
%
% .rot [real numer]: The clockwise angle expressed in radiants between the
% Y-axis and the ellipse. Default value is 0.
%
% .GK [0, 1]: The algorithm uses the Gauss-Kronrod adaptive quadrature rule
% if set to 1. It overrides opts.nc. Default is 0.
%
% .M [Positive integer]: Number of moments used by the contour integration.
% If not set, is computed at runtime.
%
% .V [Matrix]: The initial projection space used to compute the moments.
% The default value is computing the size at runtime and draw random
% unitary columns.
%
% .Vsize [Positive integer]: If V is not given, the function builds a
% random matrix with Vsize columns. If both the .V and .Vsize fields are
% given, CONTOURSOLVER ignores the latter and throws a warning. The default
% value is computed at runtime.
%
% .onlyIn [Boolean]: This parameter allow the user to retrieve the
% eigenvalues computed by CONTOURSOLVER when they are OUTSIDE the contour
% Gamma if set to 0. Theoretical resuts do not hold in this case. Default
% is 1.
%
% .Fp [funtion_handle]: The derivative of F. CONTOURSOLVER may exploit it
% to estimate the number of eigenvalues.
%
% .thres [double]: A threshold parameter used for the backward error of the
% eigenvalues. Default is 1e-11.
%
% .ref2 [0 1 2 3 4]: Choose the refinement for CONTOURSOLVER. If OPTS.ref2
% = 0, then there is no refinement; if OPTS.ref2 = 1, it is chosen
% automatically; if OPTS.ref2 = 2, then it performs a Newton refinement; if
% OPTS.ref2 = 3, then it calls itself recursively on smaller on smaller
% circles, computed with the Kmeans algorithm. If OPTS.ref2 = 4, it either
% increases the size of the projection subspace or the number of moments.
% Default is 1.
%
% .maxRefine [positive integer]: if CONTOURSOLVER does a global refinement
% by increasing the size of the projection space, this parameter sets how
% many times this refinement should be repetead. Default is 1.
%
% .NewtRefs [positive integer]: Mmaximum number of iteration for the Newton
% refinement if performed. Default is 5.
%
% .maxRecLvl [positive integer]: Maximum level of recursion if recursive
% refinement is performed. Default is 1.
%
% .method ['std', 'loewner']: This parameter chooses which method to call,
% either the standard Hankel or the Loewner interpretation. The default
% value is 'std'.
%
% .lData [struct]: If the method chosen is 'loewner', the user can set
%  the sampling points and the directions. The algorithm will then not
%  check if these parameters are good enough, therefore use this option
%  with caution. This structure must have two fields: .dir contains the 2*S
%  left and right directions; .points contains 2*S sampling points, the odd
%  indexes the left points, the even indexes the right points. By default,
%  the number S is computed during runtime, the directions are uniformly
%  independently randomly drawn,
% and the sampling points lie in a circle of radius eps^(-1/nc) centered in
% gamma.
%
% Info is a structure which contains additional information on the
% algorithm. It has the same fields , info has the following fields:
%
% .gam and .rad: the center and the radius of the target set.
%
% .w and .z: the quadrature points and the pointscwhere F is evaluated. At
% the moment, they are the nc-th roots of the unity and z = center +
% radius*w.
%
% .isInside: if the opts.onlyIn = 0, then info.isInside returns the mask of
% the eigenvalues inside the contour Gamma.
%
% .badEvs: info.badEvs is a mask such that evs(badEvs) were the eigenvalues
% whose residuals were larger than thres.
%
% .estNevs: The estimated number of eigenvalues minus poles.
%
% .err: the (potential) error message if the algorithm encounters a
% problem. It should always be an empty string.
%
% .normF: The estimation. Practically, it is computed as the max of the
% Frobenius norm of F(z)over all the integration points and the
% eigenvalues. This value is used at the denominator of the backward errors
% contained in resids.



%% EC code
expQuad = 1;
% %%

if nargout > 4
    info.err = 'Too many outputs';
    error('myfuns:contourSolver:TooManyOutputs', ...
        'At  most four outputs: eigenvalues, eigenvectors, residuals and info.');
end

if nargin > 4
    info.err = 'Too many inputs';
    error('myfuns:contourSolver:TooManyInputs', ...
        'At  most four inputs: a  matrix function F, a center, a rad, and a struct of options.');
end

%% Checking all the variables
if nargin <= 3
    [gam, rad, str, rot, n, nc, w, z, M, Vexists, V, Vsize, Fp, Fpexists, ...
        thres, maxRefine, fastUp, onlyIn, NewtRefs, ref2, nF, maxRecLvl, ...
        gk, method, lData, verbose, info] = iCheck(F, gam, rad);
else
    opts = varargin{1};
    [gam, rad, str, rot, n, nc, w, z, M, Vexists, V, Vsize, Fp, Fpexists, ...
        thres, maxRefine, fastUp, onlyIn, NewtRefs, ref2, nF, maxRecLvl, ...
        gk, method, lData, verbose, info] = iCheck(F, gam, rad, opts);
end

GK = isfield(gk, 'weights');
[lupoints, normF] = iLUfact(F, nc, z);
normF = max(normF, nF);
info.normF = normF;
%% Estimating the number of eigenvalues (eigenvalues - poles in the meromorphic case)
nevs = iContourCount(F, n, nc, rad, str, rot, z, ...
    lupoints, Fp, Fpexists, verbose);
info.estNevs = nevs;
%% Core of the algorithm
if isequal(method, 'std')
    %% Hankel method
    
    if M*Vsize < nevs
        if Vexists && verbose >=1
            % throw a warning, because the V given by is not large enough
        end
        Mnew = 1;
        Vsizenew = min(round(1.1*nevs/Mnew), n);
        Mnew = ceil(nevs/Vsizenew);
        if verbose >= 1
            warning("The given number of moments %d and size of the probing matrix %d are not large enough. They have been increased to %d and %d.", M, Vsize, Mnew, Vsizenew)
        end
        M = Mnew;
        V1 = iBuildV(n, Vsizenew-Vsize);
        V = [V V1];
        Vsize = Vsizenew;
    end
    Pspace = zeros(0,0,nc);
    Pspace = iComputeProSpace(lupoints, n, nc, V, Pspace);
    %EXPERIMENTAL CODE
    % The following is used for use to better understand how the number of
    % trapezoidal points influences the eigenvalues. I will denote with EC
    % where I need to add other pieces
    nc2 = nc*2;
    theta2 = 2*pi*(1:nc2)/nc2;
    w2 = exp(1i*theta2);
    nc4 = nc*4;
    theta4 = 2*pi*(1:nc4)/nc4;
    w4 = exp(1i*theta4);
    z2 = gam +rad*w2;
    z4 = gam + rad*w4;
    lupoints2 = iLUfact(F, nc2, z2);
    lupoints4 = iLUfact(F, nc4, z4);
    Pspace2 = zeros(0,0,nc2);
    Pspace2 = iComputeProSpace(lupoints2, n, nc2, V, Pspace2);
    Pspace4 = zeros(0,0,nc4);
    Pspace4 = iComputeProSpace(lupoints4, n, nc4, V, Pspace4);
    % End EC
    if ~GK
        Moms = iComputeMoments(Pspace, w, str, rot, nc, 0, 2*M, rad);
        %EC
        if expQuad
            Moms2 = iComputeMoments(Pspace2, w2, str, rot, nc2, 0, 2*M, rad);
            Moms4 = iComputeMoments(Pspace4, w4, str, rot, nc4, 0, 2*M, rad);
            errQuad = zeros(2,2*M);
            for j = 1:2*M
                errQuad(1,j) = norm(Moms(:,:,j) - Moms2(:,:,j),'inf')...
                    /(norm(Moms(:,:,j), 'inf'));
                errQuad(2,j) = norm(Moms2(:,:,j) - Moms4(:,:,j), 'inf')...
                    /(norm(Moms2(:,:,j), 'inf'));
            end
            %End EC
            info.errQuad = errQuad;
            info.Moms = Moms;
            info.Moms2 = Moms2;
        end
        %End EC
    else
        [Moms, ~, info.nc] = iComputeMomentsGK(Pspace, str, rot, nc, 0, 2*M, ...
            F, gam, rad, V, gk);
    end
    B0 = iBuildHank(Moms, M, 0);
    % With this syntax, svd returns a vector, not a matrix!
    sigma = svd(B0);
    tol = 1e-8; % tolerance for singular values
    mbar = find(sigma/sigma(1) > tol,1, 'last');
    if verbose >= 2
        fprintf('Performed first SVD.\n')
    end
    if nevs == 0 && mbar == Vsize*M
        if verbose >= 2
            fprintf("There are no eigenvalues inside this contour.\n")
        end
        evs = [];
        evecs = [];
        resids = [];
        return
    end
    
    VsizeMax = max(min(100,n), Vsize);
    VsizeStep = 10;
    if ~Vexists
        const = 1.2; %% a heuristic constant.
        [Pspace, Moms, V] = iUpdateV(F, n, nc, gam, rad, str, rot, w, M,...
            V, lupoints, Pspace, Moms, mbar, nevs, VsizeStep, const, ...
            Vsize, VsizeMax, tol, gk, verbose);
    else
        % We are given V, therefore we do not need to do these checkings
        if verbose >= 1
            fprintf('\nV is given, therefore we skip the estimation of the appropriate dimension.\n');
        end
    end
    M = size(Moms,3)/2;
    info.M = M;
    
    for i = 0 : maxRefine
        B0 = iBuildHank(Moms, M, 0);
        Vsize = size(B0, 2)/M;
        [Wl, Sig, Wr] = svd(B0);
        sigma = diag(Sig);
        mbar = find(sigma/sigma(1) > tol,1, 'last');
        %% Check if the projection space is good enough
        if (sigma(1) < tol || mbar == M * Vsize)
            if i ~= maxRefine
                if verbose >= 2
                    fprintf('Smallest singular value: %7.2e\n', sigma(end));
                    fprintf(['The minimum singular value is not small enough,', ...
                        ' applying refinement for the %d-th time.\n'], i+1);
                end
                %% We take A_0 as the new projection space
                if ~GK
                    Vcell = iComputeMoments(Pspace, w, str, rot, nc, 0, 1, rad);
                else
                    Vcell = iComputeMomentsGK(Pspace, str, rot, nc, 0, 1, F, gam, rad, V, gk);
                end
                V = Vcell(:,:,1);
                [~, Pspace] = iComputeProSpace(lupoints, n, nc, V, Pspace);
                continue
            else
                if verbose >= 1
                    fprintf('Smallest singular value: %7.2e\n', sigma(end));
                    warning(sprintf('%s\n%s%d%s\n', 'The minimum singular value is not small enough,', ...
                        'even after ', maxRefine, ' refinements. The results may be inaccurate'));
                end
            end
        end
        B1 = iBuildHank(Moms, M, 1);
        Wl0 = Wl(:,1:mbar);
        Sig0 = Sig(1:mbar,1:mbar);
        Wr0 = Wr(:, 1:mbar);
        Pen = (Wl0'*B1*Wr0)/Sig0;
        Wl01 = Wl0(1:n,:);
        [evecs, evs] = eig(Pen);
        info.B0 = B0;
        info.B1 = B1;
        % Eigenvectors of the NL problem
        evecs = Wl01 * evecs;
        % Eigenvalues of the NL problem
        evs = gam + rad*diag(evs);
        % Sorting them out
        [evs, evecs, resids, isInside, normF] = iCleanEigenpairs(F, ...
            evs, evecs, gam, rad, str, rot, onlyIn, normF);
        resMax = max(resids);
        if resMax <= thres*normF
            if ~onlyIn
                info.isInside = isInside;
            end
            resids = resids/normF;
            info.normF = normF;
            info.VsizeFinal = Vsize;
            return
        end
        badEvs = resids > thres*normF;
        nBadEvs = sum(badEvs);
        % Estimate the cost of the refinements
        K = 3;
        Cn = 3*n^3*nBadEvs; %3n^3b_e
        Cp = 14*n*M^3*Vsize^2 + 2*M*nc*n^2;%14nN(Ne)3/l +2pNn2;

        Cr = 2*nc*n^3/3 +14*n*M^3*Vsize^2/K^2; %2/3N n^3 +14(Ne)^3/K^2/l
        info.Cn = Cn;
        info.Cp = Cp;
        info.Cr = Cr;
        % We never want to automatically choose Cr.
        Cr = Cp+1;
        if verbose >= 1
            fprintf(['%7.0f residuals are larger than thres = %7.2e:', ...
                ' max(resids) = %7.3e \n'], nBadEvs , thres, resMax/normF);
        end
        if ref2 == 2
            Cn = -1;
        elseif ref2 == 3
            Cr = -1;
        elseif ref2 == 4
            Cp = -1;
        end
        if ref2 ~= 0
            switch min([Cp, Cn, Cr])
                
                case Cp
                    if i ~= maxRefine
                        if verbose >=1
                            fprintf(['Apply refinement by increasing the ', ...
                                'size of the Hankel matrix.\n\n']);
                        end
                        % If possible, we increase the size of the probing
                        % space
                        if Vsize < VsizeMax
                            VsizeStep = min(VsizeStep, VsizeMax-Vsize);
                            V1 = iBuildV(n, VsizeStep);
                            Pspace = iComputeProSpace(lupoints, n, nc, V1, Pspace);
                            info.V = [V, V1];
                            if ~GK
                                Moms = iComputeMoments(Pspace, w, str, rot, ...
                                    nc, 0, 2*M, rad);
                            else
                                V = [V, V1];
                                Moms = iComputeMomentsGK(Pspace, str, ...
                                    rot, nc, 0, 2*M, F, gam, rad, V, gk);
                            end
                        else % we increase the number of moments by 1
                            if M == MMax
                                if verbose >= 1
                                    warning(['The refinement stopped ',...
                                        'before max(resids) <= %4.2e, ', ...
                                        'because it reached the maximum ', ...
                                        'number of moments %d'], thres, MMax);
                                end
                                info.normF = normF;
                                resids = resids/normF;
                                return
                            else
                                Mnew = M+1;
                                if ~GK
                                    MomsNew = iComputeMoments(Pspace, w, ...
                                        str, rot, nc, 2*M, 2*MNew, rad);
                                else
                                    MomsNew = iComputeMomentsGK(Pspace, ...
                                        str, rot, nc, 2*M, 2*MNew, F, gam, ...
                                        rad, V, gk);
                                end
                                Moms = cat(3, Moms, MomsNew);
                                M=Mnew;
                                info.M = M;
                            end
                        end
                    else
                        if verbose >= 1
                            warning(['The refinement stopped before ', ...
                                'max(resids) <= %4.2e'], thres);
                        end
                        info.normF=normF;
                        resids = resids/normF;
                        return
                    end
                case Cn
                    [evsNew, evecsNew, residsNew, normF] = iNewtonRef(F, ...
                        n, evs(badEvs), evecs(:, badEvs), ...
                        resids(badEvs), thres, NewtRefs, normF, verbose);
                    info.badEvs = badEvs;
                    info.nBadEvs = nBadEvs;
                    info.normF = normF;
                    evs(badEvs) = evsNew;
                    evecs(:, badEvs) = evecsNew;
                    resids(badEvs) = residsNew;
                    % we finally divide the all the residuals by normF
                    resids = resids/normF;
                    return
                case Cr
                    if maxRecLvl == 0
                        if verbose >= 1
                            warning(["Reached max number of recursion",...
                                "and max resids higher than given threshold.\n"])
                        end
                    else
                        visual = 1;
                        points = [real(evs(badEvs)) imag(evs(badEvs))];
                        [IDX, C, D] = iChoiceOfKClusters(n, w, gam, rad,...
                            points, nBadEvs, visual);
                        K = size(D,2);
                        refinedEvs = zeros(nBadEvs,1);
                        refinedEvecs = zeros(n,nBadEvs);
                        refinedResids = zeros(nBadEvs,1);
                        radii = zeros(1,K);
                        newOpts = iBuildOpts(nc, str, rot, M, Vexists, V, ...
                            Vsize, Fp, Fpexists, thres, maxRefine, fastUp, ...
                            onlyIn, NewtRefs, ref2, maxRecLvl-1, GK, verbose);
                        flagGoodEvs = zeros(length(evs),1);
                        % We call recursively contourSolver on the K cluster. This
                        % could be parallelized in the future. Notice that we need to
                        % keep track of the eigenvalues that could be computed more than
                        % once
                        for j = 1:K
                            %cPoints(:,j) = points(IDX==j,1)+1i*points(IDX==j,2); % points in cluster j
                            distances = sqrt(D(IDX==j,j)); % euclidean distances of points in cluster j from their
                            radii(j) = max(rad*1e-5,max(distances)*1.05);
                            
                            % If in the current cluster there aren't many eigenvalues, we
                            % call the Newton refinement 
                            % New code
                            cent = C(j,1)+1i*C(j,2);
                            evsInDiskJ = iIsInsideContour(evs, cent, radii(j));
                            badEvsInDiskJ = badEvs & evsInDiskJ;
                            nBadEvsInDiskJ = sum(badEvsInDiskJ);
                            nEvsInDiskJ = sum(evsInDiskJ);
                            ncJ = 2*max(round(16*2^log10(radii(j))), 2);
                            Cn = 3*n^3*nBadEvsInDiskJ;
                            VsizeInJ = Vsize; %temporaneo
                            Cri = 2*ncJ*n^3/3 + 14*nEvsInDiskJ^3/VsizeInJ; %2/3 n^3 +14(Ne)^3/K^2/l
                            if ref2 == 3
                                Cri = -1;
                            end
                            if Cn <= Cri
                                flagGoodEvs = flagGoodEvs | badEvsInDiskJ;
                                [newEvs, newEvecs, newResids, normF] = ...
                                    iNewtonRef(F, n, ...
                                    evs(badEvsInDiskJ), evecs(:, badEvsInDiskJ), ...
                                    resids(badEvsInDiskJ), ...
                                    thres, NewtRefs, normF, verbose);
                            else
                                flagGoodEvs = flagGoodEvs | evsInDiskJ;
                                % creating the opts for the recursive alg
                                newOpts.maxRecLvl = maxRecLvl-1;
                                newOpts.nF = normF;
                                newOpts.V = [evecs(:, evsInDiskJ) rand(n,0)];
                                [newEvs, newEvecs, newResids, infor] = ...
                                  contourSolver(F, cent, radii(j), newOpts);
                                % These new residuals are scaled by the local
                                % normF, which can be different from the
                                % global one, hence we have to rescale
                                newResids = newResids*infor.normF;
                                normF = max(normF, infor.normF);
                            end
                            [refinedEvs, refinedEvecs, refinedResids] = ...
                                iCleanNewtEigenPairs(newEvs, newEvecs, ...
                                newResids, refinedEvs,refinedEvecs, ...
                                refinedResids, C, radii, j);
                        end
                        %Removing the refined eigenvalues,e'vecs, and resids
                        badEvs = badEvs | flagGoodEvs;
                        evs = evs(~badEvs);
                        evecs = evecs(:,~badEvs);
                        resids = resids(~badEvs);
                        %done
                        evs = [evs; refinedEvs];
                        evecs = [evecs refinedEvecs];
                        resids = [resids; refinedResids];
                        % We can now normalise with the new normF
                        info.normF = normF;
                        resids = resids/normF;
                        % remove eigenvalues outside the original contour
                        isInside = iIsInsideContour(evs, gam, rad, str, rot);
                        if onlyIn
                            evs = evs(isInside);
                            evecs = evecs(:, isInside);
                            resids = resids(isInside);
                        end
                        % sorting evs, evecs and resids
                        [~, idx] = sort(real(evs));
                        evs = evs(idx);
                        evecs = evecs(:, idx);
                        resids = resids(idx);
                    end
                    return
            end
        else % we do not refine and we divide the residuals by normF
            resids = resids/normF;
            return
        end
    end

elseif isequal(method, 'loewner')
    %% Loewner
    
    %% EXPERIMENTAL CODE
    % The following is used for use to better understand how the number of
    % trapezoidal points influences the eigenvalues. I will denote with EC
    % where I need to add other pieces
    nc2 = nc*2;
    theta2 = 2*pi*(1:nc2)/nc2;
    w2 = exp(1i*theta2);
    nc4 = nc*4;
    theta4 = 2*pi*(1:nc4)/nc4;
    w4 = exp(1i*theta4);
    z2 = gam +rad*w2;
    z4 = gam + rad*w4;
    lupoints2 = iLUfact(F, nc2, z2);
    lupoints4 = iLUfact(F, nc4, z4);
    %% End EC
    nLoewPoints = length(lData.points)/2;
    lDir = lData.dir(:,1:nLoewPoints); %left directions
    rDir = lData.dir(:,nLoewPoints+1:end); %right directions
    % lData.exists is true if lData was given by the user, otherwise is
    % false.
    if nLoewPoints <= nevs
        if verbose >= 1 && lData.exists
            warning(['The given number of points %d for the Loewner ', ...
                ' interpolation is not large enough. This may cause ', ...
                'inaccuracies.'], nLoewPoints)
        end
        % We do not change anything if the data is given as an input
        if ~lData.exists
            nLoewPointsNew = round(nevs*1.05)+1;
            Dir1 = iBuildV(n, 2*(nLoewPointsNew-nLoewPoints));
            auxSize = size(Dir1,2);
            lDir = [lDir, Dir1(:,1:auxSize/2)]; %left directions
            rDir = [rDir, Dir1(:,auxSize/2+1:end)]; %right directions
            nLoewPoints = nLoewPointsNew;
        end
    end
    if lData.exists
        lPoints = lData.points(1:2:end);
        rPoints = lData.points(2:2:end);
    else
        [lPoints, rPoints, waux] = iBuildLoewnerData(nLoewPoints, nc, gam,...
            rad, str, rot);
    end
    info.lPoints = lPoints;
    info.rPoints = rPoints;
    rankDeficient = 0;
    
    for j = 0:maxRefine
        [lf, rf] = iComputeLoewnerIntegral(lupoints, z, n, nc, rad,...
            str, rot, lPoints, rPoints, lDir, rDir);
        while ~rankDeficient
            % Experimental code
            if expQuad
                [lf2, rf2] = iComputeLoewnerIntegral(lupoints2, z2, n, nc2, rad, str, rot, ...
                    lPoints, rPoints, lDir, rDir);
                [lf4, rf4] = iComputeLoewnerIntegral(lupoints4, z4, n, nc4, rad, str, rot, ...
                    lPoints, rPoints, lDir, rDir);
                errQuad = zeros(2,2*nLoewPoints);
                errQuad(1,:) = vecnorm([lf-lf2, rf-rf2])./vecnorm([lf2, rf2]);
                errQuad(2,:) = vecnorm([lf2-lf4, rf2-rf4])./vecnorm([lf4, rf4]);
                info.errQuad = errQuad;
            end
            % End experimental code
            [L, Ls, F1, G1] = iBuildLoewnerMatrices(lDir, rDir, lPoints,...
                rPoints, lf, rf);
            info.L = L;
            info.Ls = Ls;
            % Check that the rank is deficient
            tol = 1e-8;
            [X, Sigma1, ~] = svd([L, Ls]);
            mbar = find(diag(Sigma1)/Sigma1(1,1) > tol,1, 'last');
            if nevs == 0 && mbar == nLoewPoints
                if verbose >= 2
                    fprintf("There are no eigenvalues inside this contour.\n")
                end
                evs = [];
                evecs = [];
                resids = [];
                return
            end
            if mbar == nLoewPoints && ~lData.exists
                % it means that we have not used enough points. We do not
                % update if we were given the data as input
                nLoewPointsNew = round((nLoewPoints+1)*1.05);
                Dir1 = iBuildV(n, 2*(nLoewPointsNew-nLoewPoints));
                auxSize = size(Dir1,2);
                lDirNew = Dir1(:,1:auxSize/2);
                rDirNew = Dir1(:,auxSize/2+1:end);
                lDir = [lDir, lDirNew]; %left directions
                rDir = [rDir, rDirNew]; %right directions
                aux = nLoewPointsNew-nLoewPoints;
                [lPointsNew, rPointsNew, waux] = iUpdateLoewnerData(aux,...
                    lPoints, rPoints, waux, gam);
                [lfNew, rfNew] = iComputeLoewnerIntegral(lupoints, z,...
                    n, nc, rad, str, rot, lPointsNew, rPointsNew, ...
                    lDirNew, rDirNew);
                lf = [lf lfNew];
                rf = [rf rfNew];
                lPoints = [lPoints lPointsNew];
                rPoints = [rPoints rPointsNew];
                nLoewPoints = nLoewPointsNew;
                info.lPoints = lPoints;
                info.rPoints = rPoints;
            else
                rankDeficient = 1;
            end
        end
        [~, ~, Y] = svd([L; Ls]);
        X0 = X(:,1:mbar);
        Y0 = Y(:,1:mbar);
        Lbar = X0'*L*Y0;
        Lsbar = X0'*Ls*Y0;
        F1bar = F1*Y0;
        G1bar = G1*X0;
        H = @(z) F1bar/(Lsbar - z*Lbar)*G1bar';
        % eigenvalues and eigenvectors
        [evecs, evs] = eig(Lsbar, Lbar);
        evs = diag(evs);
        evecs = F1bar*evecs;
        [evs, evecs, resids, isInside] = iCleanEigenpairs(F, evs, evecs,...
            gam, rad, str, rot, onlyIn, normF);
        resMax = max(resids);
        if resMax < thres*normF
            if ~onlyIn
                info.isInside = isInside;
            end
            info.VsizeFinal = Vsize;
            return
        end
        if verbose >=1
            fprintf(['Some residuals are larger than thres ',...
                '= %7.2e: max(resids) = %7.3e. \n'], thres*normF, resMax);
        end
        resids = resids/normF;
        if j ~= maxRefine
            if verbose >= 1
                fprintf(['Refinement using the eigenvectors as (some) of ', ...
                    'the right directions.\n\n']);
            end
            % We always know that #evs <= r;
            % We use the eigenvectors as newleft and right directions.
            newlength = floor(size(rDir, 2)/size(evecs,2));
            newDir = kron(ones(1,newlength), evecs);
            newDir = [newDir evecs(:,1:mod(size(rDir,2),size(evecs,2)))];
            rDir = newDir;
            lDir = newDir;            
            [lf, rf] = iComputeLoewnerIntegral(lupoints, z, n, nc,rad, ...
                str, rot, lPoints, rPoints, lDir, rDir);
            [L, Ls, F1, G1] = iBuildLoewnerMatrices(lDir, rDir, ...
                lPoints, rPoints, lf, rf);
            [X, ~, ~] = svd([L, Ls]);
        else
            if verbose >= 1
                warning(['The refinement stopped, but  max(resids) ',... 
                    '= %4.2e >= %4.2e'],...
                    resMax, thres*normF);
            end
        end
        info.L = L;
        info.Ls = Ls;
    end 
end
end



%% Auxiliary functions

function [Moms, MomsPerm] = iComputeMoments(Prospace, w, str, ...
    rot, nc, Mlow, Mup, rad, varargin)

% [Moms, MomsPerm] = COMPUTEMOMENTS(Pspace, w, nc, str, rot, nc, Mlow, Mup,
% varargin) is an auxiliary function of CONTOURSOLVER. It computes the
% moments of a n x n matrix function F(z), i.e., the integrals
%
%     A_k(V) = \int_Gamma z^k * F(z)^{-1} * V dz
%
% for k = Mlow, ..., Mup-1.
%
% INPUT
% - Pspace: It is a 'n-by-Vsize-by-nc' tensor. At each layer, it contains the
%   matrices Pspace(:,:,j) = F(z_j)^{-1}*V for j = 1, ..., nc.
%   See 'help COMPUTEPROSPACE'.
% - w: A vector of size nc with the integration points/weights.
%   They are the nc-th roots of 1, because we use the
%   trapezoidal rule. For the adaptive rule, see COMPUTEMOMENTSGK.M
% - str and rot: The stretch and rotation parameter to identify the contour
% Gamma
% - Mlow and Mup: The lowest and largest moment to compute.
% - rad: The radius of the contour Gamma.
% - varargin: if given, then it corresponds to PspaceNew, the
%   second output of COMPUTEPROSPACE.
%
% OUTPUT
% - Moms: It's a 'n-by-Vsize-by-M' tensor, Moms(:,:,k) = A_{k-1}(V) =
%   \int_Gamma z^{k-1} * F(z)^{-1} * V dz computed with the trapezoidal rule.
%   If PspaceNew is given, then Moms(:,:,k) = A_k([V V1])
% - MomsPerm: if PspaceNew is given, then MomsPerm(:,:,k) = A_{k-1}(V1).

% The differential of the integral of the ellipse, i.e. the derivative of
% w. See check.m to retrieve the formula. The 2pi is omitted due to the
% fact that it is simplified in the final integral, so in reality dw is
% "dw/2pi*i". If rot is 0, dw = w;

theta = 2*pi*(1:nc)/nc;
dw = 1i*sin(theta).*cos(rot) + 1i*str*cos(theta).*sin(rot) ...
    - sin(rot).*sin(theta) + str*cos(rot).*cos(theta);
nMoms = Mup - Mlow;

if (nargin == 9)
    ProspacePerm = varargin{1};
    [n, VsizePerm, ~] = size(ProspacePerm);
    MomsPerm = zeros(n, VsizePerm, nMoms);
    [~, Vsize, ~] = size(Prospace);
    Moms = zeros(n, Vsize+VsizePerm, nMoms);
    for k = Mlow+1:Mup
        reshaped = reshape(rad*w.^(k-1).*dw, [1,1,nc]);
        MomsPerm(:,:,k-Mlow) = sum(bsxfun(@times,ProspacePerm,reshaped),3)/nc;
        Moms(:,:,k-Mlow) = sum(bsxfun(@times,Prospace,reshaped),3)/nc;
    end
    Moms = cat(2, Moms, MomsPerm);
else
    [n,Vsize, ~] = size(Prospace);
    Moms = zeros(n,Vsize,nMoms);
    for k = Mlow+1:Mup
        reshaped = reshape(rad*w.^(k-1).*dw, [1,1,nc]);
        Moms(:,:,k-Mlow) = sum(bsxfun(@times,Prospace,reshaped),3)/nc;
    end
end
end

function [Moms, MomsPerm, nPoints] = iComputeMomentsGK(Prospace, str,...
    rot, nc, Mlow, Mup, F, gam, rad, V, gk, varargin)

% [Moms, MomsPerm, nPoints] = COMPUTEMOMENTSGK(Pspace, str, rot, nc,
% Mlow, Mup, F, gam, rad, V, gk varargin) is an auxiliary function of
% CONTOURSOLVER. It computes the moments of a n-by-n matrix function F(z),
% i.e., the integrals
%
%     A_k(V) = \int_Gamma z^k * F(z)^{-1} * V dz
%
% for k = Mlow, ..., Mup-1. It uses the Gauss-Kronrod adaptive rule.
%
% INPUT
% - Pspace: It is a 'n-by-Vsize-by-nc' tensor. At each layer, it contains
%    the matrices Pspace(:,:,j) = F(z_j)^{-1}*V for j = 1, ..., nc.
%    See 'help COMPUTEPROSPACE'.
% - w: A vector of size nc with the integration points/weights.
% - str and rot: The stretch and rotation parameters of the contour Gamma.
%    For more information, see the main help.
% - nc: The initial number of integration points.
% - Mlow and Mup: The lowest and largest moment to compute.
% - F: The matrix function given as a function-handle.
% - gam and rad: The center and the radius of the contour Gamma.
% - V: The projection subspace
% - gk: The variable that contains the points and the weights of the
%    Gauss-Kronrod quadrature rule.
% - varargin: if given, then it corresponds to PspaceNew, the
%    second output of COMPUTEPROSPACE.
%
% OUTPUT
% - Moms: It's a 'n-by-Vsize-by-M' tensor, Moms(:,:,k) = A_{k-1}(V) =
%    \int_Gamma z^{k-1} * F(z)^{-1} * V dz computed with the trapezoidal rule.
%    If PspaceNew is given, then Moms(:,:,k) = A_k([V V1])
% - MomsPerm: if PspaceNew is given, then MomsPerm(:,:,k) = A_{k-1}(V1).
% - nPoints: The final number of quadrature points used.

nMoms = Mup - Mlow;
weights = gk.weights;
weightsGauss = gk.weightsGauss;
nodes = gk.points;
threshold = 10^(-5)/2;
nPoints = 0;

if nargin == 12
    ProspacePerm = varargin{1};
    [n, VsizePerm, ~] = size(ProspacePerm);
    MomsPerm = zeros(n, VsizePerm, nMoms);
    [~, Vsize, ~] = size(Prospace);
    Moms = zeros(n, Vsize+VsizePerm, nMoms);
    aux.Prospace = Prospace;
    aux.ProspacePerm = ProspacePerm;
    aux.a = 0;
    aux.b = 1;
    integralPoints = {aux};
    finished = 0;
    j = 1;
    q = 0;
    qperm = 0;
    while ~finished
        a = integralPoints{j}.a;
        b = integralPoints{j}.b;
        Prospace = integralPoints{j}.Prospace;
        ProspacePerm = integralPoints{j}.ProspacePerm;
        halfh = (b-a)/2;
        midpt = (b+a)/2;
        weightsInt = weights*halfh;
        weightsGaussInt = weightsGauss*halfh;
        ww = 2*pi*(nodes*halfh + midpt);
        wwEllipse = cos(ww)*cos(rot) - str*sin(ww)*sin(rot) + ...
            1i*sin(rot)*cos(ww) + 1i*str*cos(rot)*sin(ww);
        dw = 1i*sin(ww).*cos(rot) + 1i*str*cos(ww).*sin(rot) - ...
            sin(rot).*sin(ww) + str*cos(rot).*cos(ww);
        reshaped = reshape(rad*wwEllipse.^(Mup-1).*dw.*weightsInt, [1,1,nc]);
        reshapedGauss = reshape(rad*wwEllipse(2:2:end).^(Mup-1).* ...
            dw(2:2:end).*weightsGaussInt, [1,1,floor(nc/2)]);
        qtemp = sum(bsxfun(@times,Prospace,reshaped),3);
        qtempPerm = sum(bsxfun(@times,ProspacePerm,reshaped),3);
        qtempGaussPerm = sum(bsxfun(@times,ProspacePerm(:,:,2:2:end), ...
            reshapedGauss),3);
        errtemp = norm(qtempPerm-qtempGaussPerm);
        if (errtemp < threshold || halfh <= 2^(-6))
            qperm = qperm + qtempPerm;
            q = q + qtemp;
            errtot = errtot + ertemp;
            integralPoints{j}.used = 1;
        else
            aux1.a = a; aux1.b = midpt;
            aux2.a = midpt; aux2.b = b;
            newPoints1 = 2*pi*(nodes*halfh/2 +(midpt+a)/2);
            newPoints2 = 2*pi*(nodes*halfh/2 +(midpt+b)/2);
            newPoints1Ellipse = cos(newPoints1)*cos(rot) - ...
                str*sin(newPoints1)*sin(rot) + 1i*sin(rot)*cos(newPoints1) ...
                + 1i*str*cos(rot)*sin(newPoints1);
            newPoints2Ellipse = cos(newPoints2)*cos(rot) - ...
                str*sin(newPoints2)*sin(rot) + 1i*sin(rot)*cos(newPoints2) ...
                + 1i*str*cos(rot)*sin(newPoints2);
            zz1 = rad*newPoints1Ellipse + gam;
            zz2 = rad*newPoints2Ellipse + gam;
            lupoints = iLUfact(F,nc,zz1);
            aux1.Prospace = iComputeProSpace(lupoints, n, nc, V, []);
            [~,aux1.ProspacePerm] = iComputeProSpace(lupoints, n ,nc, V1, aux1.Prospace);
            lupoints = iLUfact(F,nc,zz2);
            aux2.Prospace = iComputeProSpace(lupoints, n, nc, V, []);
            [~,aux2.ProspacePerm] = iComputeProSpace(lupoints, n ,nc, V1, aux2.Prospace);
            integralPoints(end+1:end+2) = {aux1 aux2};
            integralPoints{j}.used = 0;
        end
        if j == length(integralPoints)
            finished = 1;
        end
        j = j+1;
    end
    MomsPerm(:,:,nMoms) = qperm;
    Moms(:,:,nMoms) = q;
    
    for k = Mlow+1:Mup-1
        q = 0;
        qperm = 0;
        for j = 1:length(integralPoints)
            if integralPoints{j}.used
                a = integralPoints{j}.a;
                b = integralPoints{j}.b;
                Prospace = integralPoints.Prospace;
                ProspacePerm = integralPoints.ProspacePerm;Un os
                midpt = (b+a)/2;
                halfh = (b-a)/2;
                weightsInt = weights*halfh;
                ww = 2*pi*(nodes*halfh + midpt);
                wwEllipse = cos(ww)*cos(rot) - str*sin(ww)*sin(rot) + ...
                    1i*sin(rot)*cos(ww) + 1i*str*cos(rot)*sin(ww);
                dw = 1i*sin(ww).*cos(rot) + 1i*str*cos(ww).*sin(rot) ...
                    - sin(rot).*sin(ww) + str*cos(rot).*cos(ww);
                reshaped = reshape(rad*wwEllipse.^(k-1).*dw.*weightsInt, [1,1,nc]);
                qtempperm = sum(bsxfun(@times,ProspacePerm,reshaped),3);
                qtemp = sum(bsxfun(@times,Prospace,reshaped),3);
                q = q + qtemp;
                qperm = qperm + qtempperm;
                nPoints = nPoints+length(nodes);
            end
        end
        MomsPerm(:,:,k-Mlow) = qperm;
        Moms(:,:,k-Mlow) = q;
    end
    Moms = cat(2, Moms, MomsPerm);
    
else
    [n,Vsize, ~] = size(Prospace);
    Moms = zeros(n,Vsize,nMoms);
    MomsGauss = Moms;
    
    % k = M
    % Initialize main loop
    
    % IntegralPoints is a cell containing structures  with the following
    % fields:
    % - .a: left
    % - .b: right
    % - .Prospace: the matrix F(z)^(-1)V
    aux.Prospace = Prospace;
    aux.a = 0;
    aux.b = 1;
    integralPoints = {aux};
    finished = 0;
    j = 1;
    % This will contain the numerical value of the quadrature
    q = 0;
    errtot = 0;
    %First cycle to compute the zeroth Moment
    while ~finished
        a = integralPoints{j}.a;
        b = integralPoints{j}.b;
        Prospace = integralPoints{j}.Prospace;
        % half the length and middle point of the interval [a,b]
        halfh = (b-a)/2;
        midpt = (b+a)/2;
        % Getting the right weights
        weightsInt = weights*halfh;
        weightsGaussInt = weightsGauss*halfh;
        ww = 2*pi*(nodes*halfh + midpt);
        wwEllipse = cos(ww)*cos(rot) - str*sin(ww)*sin(rot) + ...
            1i*sin(rot)*cos(ww) + 1i*str*cos(rot)*sin(ww);
        dw = 1i*sin(ww).*cos(rot) + 1i*str*cos(ww).*sin(rot) ...
            - sin(rot).*sin(ww) + str*cos(rot).*cos(ww);
        % Reshaping for using bsxfun
        reshaped = reshape(rad*wwEllipse.^(Mup-1).*dw.*weightsInt, [1,1,nc]);
        reshapedGauss = reshape(rad*wwEllipse(2:2:end).^(Mup-1) ...
            .*dw(2:2:end).*weightsGaussInt, [1,1,floor(nc/2)]);
        
        % Temporary numerical value of the integral
        qtemp = sum(bsxfun(@times,Prospace,reshaped),3);
        qtempGauss = sum(bsxfun(@times,Prospace(:,:,2:2:end),reshapedGauss),3);
        % Temporary error
        errtemp = norm(qtemp-qtempGauss);
        % If the approximation is good enough, accept it and state that we
        % have used intervalPoints{j} was used
        if (errtemp < threshold*2*halfh || halfh <= 2^(-8))
            q = q + qtemp;
            % This will be the final error
            errtot = errtot + errtemp;
            integralPoints{j}.used = 1;
        else % We halve the two subinterval aux1 and aux2
            aux1.a = a; aux1.b = midpt;
            aux2.a = midpt; aux2.b = b;
            
            % We need to compute the points to compute the new Prospace
            newPoints1 = 2*pi*(nodes*halfh/2 +(midpt+a)/2);
            newPoints2 = 2*pi*(nodes*halfh/2 +(midpt+b)/2);
            newPoints1Ellipse = cos(newPoints1)*cos(rot) - ...
                str*sin(newPoints1)*sin(rot) + 1i*sin(rot)*cos(newPoints1) ...
                + 1i*str*cos(rot)*sin(newPoints1);
            newPoints2Ellipse = cos(newPoints2)*cos(rot) - ...
                str*sin(newPoints2)*sin(rot) + 1i*sin(rot)*cos(newPoints2) ...
                + 1i*str*cos(rot)*sin(newPoints2);
            zz1 = rad*newPoints1Ellipse + gam;
            zz2 = rad*newPoints2Ellipse + gam;
            lupoints = iLUfact(F,nc,zz1);
            aux1.Prospace = iComputeProSpace(lupoints, n, nc, V, []);
            lupoints = iLUfact(F,nc,zz2);
            aux2.Prospace = iComputeProSpace(lupoints, n, nc, V, []);
            integralPoints(end+1:end+2) = {aux1 aux2};
            % We set that the interval j is not used
            integralPoints{j}.used = 0;
        end
        % If we have arrived at the end, we have finished.
        if j == length(integralPoints)
            finished = 1;
        end
        j = j+1;
    end
    Moms(:,:,nMoms) = q;
    
    % Cycle to compute the other moments. These cycles are easier
    for k = Mlow+1:Mup-1
        q = 0;
        for j = 1:length(integralPoints)
            if integralPoints{j}.used
                a = integralPoints{j}.a;
                b = integralPoints{j}.b;
                Prospace = integralPoints{j}.Prospace;
                midpt = (b+a)/2;
                halfh = (b-a)/2;
                weightsInt = weights*halfh;
                ww = 2*pi*(nodes*halfh + midpt);
                wwEllipse = cos(ww)*cos(rot) - str*sin(ww)*sin(rot) + ...
                    1i*sin(rot)*cos(ww) + 1i*str*cos(rot)*sin(ww);
                dw = 1i*sin(ww).*cos(rot) + 1i*str*cos(ww).*sin(rot) ...
                    - sin(rot).*sin(ww) + str*cos(rot).*cos(ww);
                reshaped = reshape(rad*wwEllipse.^(k-1).*dw.*weightsInt, [1,1,nc]);
                qtemp = sum(bsxfun(@times,Prospace,reshaped),3);
                q = q + qtemp;
                if k == 1
                    nPoints = nPoints+length(nodes);
                end
            end
        end
        Moms(:,:,k-Mlow) = q;
    end
    %This is just to avoid errors. It is not used
    MomsPerm = 0;
    
end
end

function [ProSpace, ProSpacePerm] = iComputeProSpace(lupoints, n, nc, V, ProSpace)

% [Pspace, PspaceNew] = ICOMPUTEPROSPACE(lupoints, n, nc, V, PspaceIn)
% is an auxiliary function of CONTOURSOLVER. It is the first step to compute the
% projection space of the main algorithm.
%
% INPUT
% - lupoints: is a cell of length nc containing the LU
% factorization of F(z_j) at the integration points z_j, for j = 1,
% ..., nc. This cell is built by LUFACT.
% - n: size of F(z).
% - nc: The number of integration points.
% - V: The projection space.
% - PspaceIn: is a structure of length nc with the previous output of
% this function.
%
% OUTPUT
% - Pspace: is a structure of length nc containing the matrices
% [PspaceIn{j} PspaceNew{j}] for j = 1, ..., nc.
% - PspaceNew: is a structure of length nc with the matrices
% PspaceNew{j} = F(z_j)^{-1} * V, for j = 1, ..., nc.


l = size(V, 2);
ProSpacePerm = zeros(n, l, nc);
for j = 1:nc
    ProSpacePerm(:,:,j) = lupoints{j}\V;
end
if isempty(ProSpace)
    ProSpace = ProSpacePerm;
else
    ProSpace = cat(2, ProSpace, ProSpacePerm);
end

end

function [lf, rf] = iComputeLoewnerIntegral(lupoints, z, n, nc, ...
    rad, str, rot, lPoints, rPoints, lDir, rDir) 

% [lf, rf] = ICOMPUTELOEWNERINTEGRAL(lupoints, w, z, n, nc, rad, str, rot,
% lPoints, rPoints, lDir, rDir) is an auxiliary function of CONTOURSOLVER.
% It computes the left and right Loewner data (lf, rf) in the Loewner
% interpretation of Beyn's algorithm. The input are the same of
% ICOMPUTEMOMENTS, except for lPoints, rPoints and lDir, and rDir, which
% are the interpolation points and directions.



% The dw is in reality dw/(2*pi*i)
theta = 2*pi*(1:nc)/nc;
dw = 1i*sin(theta).*cos(rot) + 1i*str*cos(theta).*sin(rot) ...
    - sin(rot).*sin(theta) + str*cos(rot).*cos(theta);

r = length(rPoints);
lf = zeros(n,r);
rf = lf;
for j = 1:r
    rproSpace = iComputeProspaceLoewner(lupoints, n, nc, rDir(:,j), 0);
    lproSpace = iComputeProspaceLoewner(lupoints, n, nc, lDir(:,j), 1);
    rreshaped = reshape(rad*dw./(rPoints(j)-z), [1,1,nc]);
    lreshaped = reshape(conj(rad*dw./(lPoints(j)-z)), [1,1,nc]);
    rf(:,j) = sum(bsxfun(@times,rproSpace,rreshaped),3)/nc;
    lf(:,j) = sum(bsxfun(@times,lproSpace,lreshaped),3)/nc;
end

end

function proSpace = iComputeProspaceLoewner(lupoints, n, nc, V, transposeFlag)

% [lf, rf] = ICOMPUTEPROSPACELOEWNER(lupoints, n, nc, rad, V,
% transposeFlag) is an auxiliary function of CONTOURSOLVER. It corresponds
% to iCOMPUTEPROSPACE, but in the Loewner interpretation, and in fact the
% inputs and output are the same. There is an additional input,
% transposeFlag, since the Loewner integral is slightly different when
% computing the right or left interpolation data and requires a
% transposition.

l = size(V, 2);
proSpace = zeros(n, l, nc);
if ~transposeFlag
    for j = 1:nc
        proSpace(:,:,j) = lupoints{j}\V;
    end
else
    for j = 1:nc
        proSpace(:,:,j) = (V'/lupoints{j})';
    end
end


end

function [lupoints, normF] = iLUfact(F, nc, z)

% [lupoints, normF] = LUFACT(F, nc, z) is an auxiliary function for
% CONTOURSOLVER.
%
% INPUT
% - F: is a function_handle.
% - nc: is the number of points where the user evaluates the
% function.
% - z: is a vector of length
% nc containing the points where the user evaluates F.
%
% OUTPUT
% - lupoints: is a cell of nc structures containing the
% decomposition of F(z(k)). We used the command decompose, which is
% available since MATLAB R2017b.
% - normF: is max_k ||F(z(k)||, where ||.|| is the Frobenius norm.

lupoints = cell(1, nc);
normF = 0;

for k = 1:nc
    Fzk = F(z(k));
    normF = max(normF, norm(Fzk, 'fro'));
    lupoints{k} = decomposition(Fzk);
end

end

function nevs = iContourCount(F, n, nc, rad, str, rot, z, ...
    lupoints, Fp, Fpexists, verbose)


% nevs  = ICONTOURCOUNT(F, n, nc, gam, rad, str, rot, w, z, lupoints, Fp,
% Fpexists, verbose) is an auxiliary function of
% CONTOURSOLVER. It  computes an approximation of the number
% of eigenvalues minus the number of poles of a meromorphic
% matrix function F(z) inside a circular region Gamma. if F(z)
% is holomorphic, then CONTOURCOUNT returns an approximation
% of the number of eigenvalues.
%
% The script is based on the formula
%
%   nevs = \frac{1}{2pi*i}*\int_Gamma trace(F(z)^{-1}F(z)') dz
%
% The trace operator is stochastically estimated.
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
% OUTPUT - nevs: the number of eigenvalue minus poles of the
% function F(z).


% The differential of the integral of the ellipse, i.e. the derivative of
% w. See check.m to retrieve the formula. The 2pi is omitted due to the
% fact that it is simplified in the final integral, so in reality dw is
% "dw/2pi*i". If rot is 0, dw = w;
theta = 2*pi*(1:nc)/nc;
dw = 1i*sin(theta).*cos(rot) + 1i*str*cos(theta).*sin(rot) ...
    - sin(rot).*sin(theta) + str*cos(rot).*cos(theta);

% number of random vectorsin stochastic estimation
L = max(50, floor(n/10));
% scaling the direction of the derivative
delta = 1e-4;
% If the matrix is small, then we compute the trace exactly by the
% formula trace(A) = sum_{i=1}^n(e_i^T*A*e_i). Otherwise, we use Rademacher
% vectors
if n <= 100
    Stovecs = eye(n);
    exact = 1;
else
    Stovecs = sign(rand(n,L) - .5);
    exact = 0;
end
% This is the final result
nr = 0;
% We consider the cases of having the derivative or not
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
        nr = nr + rad*dw(k)*estTrace;
    end
else
    for k = 1:nc
        estTrace = 0;
        for j = 1:size(Stovecs,2)
            aux =  Stovecs(:,j)'/lupoints{k};
            estTrace = estTrace + aux * (Fp(z(k))*Stovecs(:,j));
        end
        nr = nr + rad*dw(k)*estTrace;
    end
end
if ~exact
    nr = nr/L;
end
nevs = round(real(nr/nc));
if (verbose >= 1)
    fprintf(['Stochastic estimation of #(eigenvalues - '...
        'poles) with %d integration points and ' ...
        '%d vectors: %d\n'], ...
        nc, L, nevs);
end
end

function [L, Ls, F1, G1] = iBuildLoewnerMatrices(lDir, rDir, lPoints, ...
    rPoints, lf, rf)
% This simple auxiliary function build the Loewner matrices thanks to the
% Loewner data.
F1 = rf;
G1 = lf;

% Building Loewner matrices
CC = 1./(lPoints.'-rPoints); % Cauchy matrix
L = (lf'*rDir - lDir'*rf).*CC;
Ls = (diag(lPoints)*lf'*rDir - lDir'*rf*diag(rPoints)).*CC;

end

function [lPoints, rPoints, waux] = iBuildLoewnerData(nLoewPoints, nc, ...
    gam, rad, str, rot)

% The default behaviour is returning the points on a arc of degree
% [-pi/4,pi/4] of  a bigger circle/ellipse than the contour one and
% putting the left and right points interlaced. The radius of the circle is
% chosen such that the eigenvalues outside are filtered out.

theta = 2*pi*(0:16*nLoewPoints-1)/(16*nLoewPoints)-pi/4;
w = cos(theta)*cos(rot) - str*sin(theta)*sin(rot) + 1i*sin(rot)*cos(theta) ...
    + 1i*str*cos(rot)*sin(theta);
% This is used if/when we need to add other interpolation points, so we return it
waux = gam + eps^(-1/nc)*rad*w(2:2:end); 
interPoints = waux(1:2*nLoewPoints);
lPoints = interPoints(1:2:end);
rPoints = interPoints(2:2:end);
end

function [lPointsNew, rPointsNew, waux] = iUpdateLoewnerData(...
    nLoewPointsNew, lPoints, rPoints, waux, gam)

% [lPointsNew, rPointsNew, waux] = iUpdateLoewnerData(nLoewPointsNew,
% lPoints, rPoints, waux, gam, rad) is an auxiliary function of
% CONTOURSOLVER, which is used in the Loewner interpretation. It is called
% when the chosen number of interpolation points is not large enough to
% compute the eigenvalues accurately. If that is the case, then points are
% added on the same circle at the same distances. If there is not enough
% space, then we add them on a vertical line on the right of the circle.

if ~isempty(waux)
    startIndex = find(rPoints(end)==waux);
    nPointsAvail = length(waux)-startIndex;
    if 2*nLoewPointsNew <= nPointsAvail
        % If we can still add points in the circle created in
        % iBuildLoewnerData, we do that
        lPointsNew = waux(startIndex+1:2:startIndex+2*nLoewPointsNew);
        rPointsNew = waux(startIndex+2:2:startIndex+2*nLoewPointsNew);
    else
        % If not, we complete the circle, then we add points on vertical line
        % indefinetely
        nPointsStillAdd = 2*nLoewPointsNew-nPointsAvail;
        if nPointsAvail % if startIndex is at the end, the next two lines do not work
            lPointsNew = waux(startIndex+1:2:end);
            rPointsNew = waux(startIndex+2:2:end);
        else
            lPointsNew = [];
            rPointsNew = [];
        end
        radius = abs(lPoints(1)-gam)*1.1;
        dist = abs(rPoints(1)-lPoints(1));
        pointsLine = gam + radius + [0:nPointsStillAdd-1]*1i*dist;
        lPointsNew = [lPointsNew pointsLine(1:2:end)];
        rPointsNew = [rPointsNew pointsLine(2:2:end)];
        waux = [];
    end
    
else
    % We have already started adding points on the vertical line
    dist = abs(rPoints(end)-lPoints(end));
    PointsNew = rPoints(end)+[1:2*nLoewPointsNew]*dist*1i;
    lPointsNew = PointsNew(1:2:end);
    rPointsNew = PointsNew(2:2:end);
end

end

function [evs, evecs, resids, isInside, normF] = ...
    iCleanEigenpairs(F, evs, evecs, gam, rad, str, rot, onlyIn, normF)
% [evs, evecs, resids, isInside, nEvsIn, normF] = iCleanEigenpairs(F, evs,
% evecs, gam, rad, str, rot, onlyIn, normF) is an auxiliary function of
% CONTOURSOLVER. It takes care of cleaning the eigenpairs return by the
% projected linear problem and transform them in the eigenpairs of the
% original problem F(lambda)v = 0. It also computes a lower bound of
% ||F||_Omega = sup_{z\in Omega} ||F(z)|| and filters the eigenvalues that
% are outside the target region when onlyIn is set to 1.

% Sort eigenvalues out
[~, idx] = sort(real(evs));
evs = evs(idx);
evecs = evecs(:, idx);
isInside = iIsInsideContour(evs, gam, rad, str, rot);
if onlyIn
    evs = evs(isInside);
    evecs = evecs(:, isInside);
end
% Normalize eigenvectors
evecs = evecs./vecnorm(evecs);
% Compute the residues
resids = zeros(size(evs));
for j = 1:length(evs)
    if abs(evs(j)) < Inf % for the corner case
        normF = max(normF, norm(F(evs(j)), 'fro'));
        resids(j) =  norm(F(evs(j))*evecs(:, j));
    else
        resids(j) = Inf;
    end
end
% The residuals, which are an upper bound of the backward errors, are
% normalised in the main code, i.e., are divided later by normF
end

function [Pspace, Moms, V] = iUpdateV(F, n, nc, gam, rad, str, rot, w, ...
    M, V,lupoints, Pspace, Moms, mbar, nevs, VsizeStep, const, Vsize, ...
    VsizeMax, tol, gk, verbose)

% [Pspace, Moms, V] = iUpdateV(F, n, nc, gam, rad, str, rot, w, ...
%    z, M, V,lupoints, Pspace, Moms, mbar, nevs, VsizeStep, const, Vsize,...
%    VsizeMax, tol, gk, verbose)
%
% is an auxiliary function CONTOURSOLVER.  It  updates the size of the
% projection space V and the number of moments M. We privilege increasing
% the size of the projection space V, when it is not possible anymore, we
% increase the number of moments M. If fastUp = 0, then UPDATEV uses an SVD
% to compute the rank of the matrix B_0(V) until it has some
% rank-deficiency; if fastUp = 1, it  uses a fast update instead. This
% latter update is left for future works.
%
% See also CONTOURSOLVER.

GK = isfield(gk, 'weights');

% We need to check whether the current B0 has full rank.
if mbar == M * Vsize
    isRankFull = 1;
else
    isRankFull = 0;
end

% If Vsize is greater than VsizeMax, it means that we need more moments
% for sure. Therefore we increase the number of moments and set VsizeNew =
% VsizeMax.
if isRankFull && (Vsize >= VsizeMax)
    MNew = ceil(const*nevs/VsizeMax); % melius abundare quam deficiere
    isRankFull = 1;
end

if (verbose >= 2)
    fprintf(['Computing a right size for V with the SVD', ...
        ' method\n'])
end

VsizeNew = min(Vsize + VsizeStep, VsizeMax);
while (isRankFull && (VsizeNew <= VsizeMax))
    VsizeAdd = VsizeNew - Vsize;
    V1 = iBuildV(n, VsizeAdd);
    if verbose >= 2
        fprintf(['Increasing the size of the projection', ...
            ' space from %d to %d\n'], Vsize, VsizeNew)
    end
    Pspace = iComputeProSpace(lupoints, n, nc, V1, Pspace);
    if ~GK
        Moms = iComputeMoments(Pspace, w, str, rot, nc, 0, 2*M, rad);
    else
        V = [V V1];
        Moms = iComputeMomentsGK(Pspace, str, rot, nc, 0, 2*M, F, gam, ...
            rad, V, gk);
    end
    B0 = iBuildHank(Moms, M, 0);
    % Check the rank with an svd
    sigma = svd(B0);
    mbar = find(sigma/sigma(1) > tol,1, 'last');
    if mbar == M * VsizeNew
        Vsize = VsizeNew;
        VsizeNew = Vsize + VsizeStep;
        if ((Vsize == VsizeMax) || (VsizeNew > VsizeMax))
            if verbose >=1
                warning(['The block Hankel matrix B0V is still ', ...
                    'full-rank, even if, the size of ', ...
                    'V is  >= than VsizeMax = %d'], VsizeMax);
            end
            MNew = M+1;
        end
    else
        isRankFull = 0;
        Vsize = VsizeNew;
    end
end

MMax = 20;
while isRankFull && (M <= MMax)
    if verbose
        fprintf(['Increasing the number of moments', ...
            ' from %d to %d\n'], M, MNew)
    end
    if ~GK
        MomsNew = iComputeMoments(Pspace, w, str, rot, nc, 2*M, 2*MNew, ...
            rad);
    else
        MomsNew = iComputeMomentsGK(Pspace, str, rot, nc, 2*M, 2*MNew, ...
            F, gam, rad, V, gk);
    end
    Moms = cat(3, Moms, MomsNew);
    B0 = iBuildHank(Moms, MNew, 0);
    % Check the rank with an svd
    sigma = svd(B0);
    mbar = find(sigma/sigma(1) > tol,1, 'last');
    if mbar == MNew * min(Vsize, n)
        M = MNew;
        MNew = M+1;
        if (MNew > MMax)
            if verbose
                warning(['The block Hankel matrix B0V is still ', ...
                    'full-rank, even if the number of moments ', ...
                    'M is  >= than MMax = %d'], MMax);
            end
        end
    else
        isRankFull = 0;
    end
end
end

function V = iBuildV(n, Vsize)

% V = BUILDV(n, Vsize) builds a n x Vsize matrix with uniformly
% distributed random columns.

V = rand(n,Vsize) - 0.5;
V = V./vecnorm(V);
end

function B =  iBuildHank(S, M, s)

% B = BUILDHANK(S, M, s) builds a square block-Hankel matrix. It is
% an auxiliary function of CONTOURSOLVER.
%
% INPUT
% - S: is a a tensor of appropriate size.
% - M: is the number of blocks of the matrix B.
% - s: If s = 0, we get B0, if s = 1, we get B1.

% OUTPUT
% - B: A block Hankel matrix such that B_{i+j-2+s} = S(:,:,i+j-1+s)
% for i,j = 1,..., M.


[n1,n2,n3] = size(S);
S1 = mat2cell(S, n1, n2, ones(n3,1));
b = hankel(1+s:M+s, M+s:2*M-1+s);
B = cell2mat(S1(b));
end

function flag = iIsInsideContour(z, center, radius, varargin)
% This short scripts checks wheter the points z are inside the circle of
% center "center" and radius "radius". The varargin is used in the
% elliptical case, where "str" is the ratio vertical axis/horizontal axis
% and "rot" is the clockwise rotation from the Y-axis.

str = 1;
rot = 0;
if nargin >= 4
    str = varargin{1};
    if nargin >= 5
        rot = varargin{2};
    end
end
w = (z - center)/radius;
flag = (real(w)*cos(rot) + imag(w)*sin(rot)).^2 + ...
    (real(w)*sin(rot) - imag(w)*cos(rot)).^2/str^2 <= 1;
end

%% Refinements auxiliary functions

function [cleanedEvs, cleanedEvecs, cleanedResids] = ...
    iCleanNewtEigenPairs(newEvs, newEvecs, newResids, refinedEvs,...
    refinedEvecs, refinedResids, C, radii, j)
%This is a helper function for contourSolver. It should not be used on its
%own. [cleanedEvs, cleanedEvecs, cleanedResids] = ICLEANNEWTEIGENPAIRS(newEvs,
%newEvecs, newResids, refinedEvs, refinedEvecs, refinedResids, C, radii, j)
%takes a set of eigenpairs with their residuals (newEvs, newEvecs and
%newResiduals) and add them to a previous set (refinedEvs, refinedEvecs,
%refinedEvs) if they were not previously computed.
% We assume there are j circles such that refinedEvs are inside the first
% j-1 circles centered in C(1:j-1) with radii(1:j-1), therefore we use this
% information to check that each newEvs is inside one of the previous
% circles. If it's not, we simply add it, otherwise we look for it n the
% set refinedEvs. Theory tells us it should be there (if not, we add it):
% we check whether the new residual is better than the previous one. If it
% is, we swap the new one with the old one.

nn = 1;
thresh = 1e-4;
Ne = length(newEvs);
cleanedEvs = zeros(Ne,1);
cleanedEvecs = zeros(size(newEvecs));
cleanedResids = zeros(Ne,1);
if j == 1
    % It is the first time we call this function, which means newEvs is the
    % first sets of refined eigenvalues and there are no previously
    % computed ones...
    cleanedEvs(1:Ne) = newEvs;
    cleanedEvecs(:,1:Ne) = newEvecs;
    cleanedResids(1:Ne) = newResids;
else
    tempEvs = cleanedEvs;
    tempEvecs = cleanedEvecs;
    tempResids = cleanedResids;
    %... there are previously computed eigenvalues. We need to check that
    %the new eigenvalues weren't computed in previous clusters. Therefore
    %we cycle through all the newEvs
    for k = 1 : Ne
        ev = newEvs(k);
        evec = newEvecs(:,k);
        resid = newResids(k);
        alreadyComp = 0;
        % is ev in a prevoius cluster? If so, it was (almost surely)
        % already computed.
        for  l = 1 : j-1
            if iIsInsideContour(ev, C(l,1)+1i*C(l,2), radii(l))
                alreadyComp = 1;
                break
            end
        end
        
        if alreadyComp
            % If the eval is inside any of the previous cluster, then we
            % look for it. We also consider the case where the eval falls
            % inside another contour, but it wasn't previously computed.
            % This should not happen, but it's better to be sure. We do
            % this with the variable foundIt.
            foundIt = 0;
            l = 1;
            % We look for the previously computed eval. If we find it, we
            % substitute it if the new residual is better than the previous
            % one.
            while ((l <= length(refinedEvs)) && ~foundIt)
                % We check if it the same eigenvalue
                sameEv = abs(ev - refinedEvs(l)) <= abs(ev)*thresh;
                % We check if it is the same eigenvector. Given that they
                % may differ by a complex number of norm 1, we need to
                % normalize. The simple idea of
                % refinedEvecs(:,l)*evec(1)/refinedEvecs(1,l) can fail (it
                % HAD!) when the first element is close to zero. So we
                % divide by the maximum element.
                [~, maxInd] = max(abs(evec));
                sameEvec = norm(evec - refinedEvecs(:,l) * ...
                    evec(maxInd)/refinedEvecs(maxInd,l)) <= thresh;
                if sameEv && sameEvec
                    % We keep the one with the better residual
                    if resid < refinedResids(l)
                        refinedEvs(l) = ev;
                        refinedEvecs(:,l) = evec;
                        refinedResids(l) = resid;
                    end
                    % else, we do nothing
                    foundIt = 1;
                end
                l = l+1;
            end
            
            % The following code should always fail (see comments above)
            if ~foundIt
                tempEvs(nn) = ev;
                tempEvecs(:,nn) = evec;
                tempResids(nn) = resid;
                nn = nn+1;
            end
        else
            tempEvs(nn) = ev;
            tempEvecs(:,nn) = evec;
            tempResids(nn) = resid;
            nn = nn+1;
        end       
    end    
    % Output
    cleanedEvs = [refinedEvs; tempEvs(1:nn-1)];
    cleanedEvecs = [refinedEvecs tempEvecs(:,1:nn-1)];
    cleanedResids = [refinedResids; tempResids(1:nn-1)];
end



end

function opts = iBuildOpts(nc, str, rot, M, Vexists, V, Vsize, Fp, Fpexists, ...
    thres, maxRefine, fastUp,  onlyIn, NewtRefs, ref2, maxRecLvl, ...
    GK, verbose)
%This function is used internally and builds the structure opts starting
%from specific inputs. It is called when the algorithms uses the recursive
%refinement.
opts.nc = nc;
opts.str = str;
opts.rot = rot;
opts.M = M;
opts.Vexists = Vexists;
opts.V = V;
opts.Vsize = Vsize;
opts.Fp = Fp;
opts.Fpexists = Fpexists;
opts.thres = thres;
opts.maxRefine = maxRefine;
opts.fastUp = fastUp;
opts.onlyIn = onlyIn;
opts.newtRefs = NewtRefs;
opts.ref2 = ref2;
opts.maxRecLvl = maxRecLvl;
opts.GK = GK;
opts.verbose = verbose;

end

function [evs, evecs, resids, normF] = iNewtonRef(F, n, evs, evecs, ...
    resids, thres, NewtRefs, normF, verbose)

%[evs, evecs, resids, normF] = iNewtonRef(F, Fp, Fpexists, n, evs, evecs,
%resids, thres, NewtRefs, normF, verbose) performs the Newton refinement on
%the eigenpairs of F(z). The only particularity is that we use a robust
%deflation appeared for the first time on "C. Effenberger. “Robust
%successive computation of eigenpairs for nonlinear  eigenvalue problems.”
%In: SIAM J. Matrix Anal. Appl. 34.3 (2013), pp. 1231–1256.". 
% In our experience, we never found an example where the approximation
% coming from the contour solver was not good enough to fall in its
% appropriate eigenvalue.

Nevs = length(evs);
u = zeros(n, 1);

delta = 1e-5; % derivative approx

% A structure that contains the minimal pairs (X, Lambda)
pairs = struct();
pairs.lambdas = cell(1,Nevs);
pairs.Xs = cell(1, Nevs);
pairs.minIndeces = ones(1,Nevs);
pairs.length = 0;
tol1 = 1e-8; %minimal distance between two eigenvalues for them
%to be equal.

for j = 1 : Nevs
    evec = evecs(:, j);
    ev = evs(j);
    [~, argmax] = max(evec);
    u(argmax) = 1;
    [ev, evec, resid, normF] = iRefinement(F, ev, evec, u, thres, ...
        delta, NewtRefs, normF, ...
        verbose);
    
    % Check if we have found a new eigenpair
    [NewEigenpair, epairIndex] = iCheckNewEig(ev, evec, pairs, tol1);
    % If it is new, we update the minimal invariant pairs.
    if NewEigenpair
        pairs.length = pairs.length + 1;
        len = pairs.length;
        pairs.lambdas{len} = ev;
        pairs.Xs{len} = evec;
    end
    
    % If it is not new, we apply the robust deflation proposed by
    % Effemberger
    while ~NewEigenpair
        %% Build A bigger minimal pair
        fprintf('This eigenpair was already computed. Applying deflation\n');
        Lambda = pairs.lambdas{epairIndex};
        X = pairs.Xs{epairIndex};
        minIndex = pairs.minIndeces(epairIndex);
        [A,B,U] = iBuildAuxFun(pairs, epairIndex, F);
        Ftil = @(z) [F(z) U(z); A(z) B(z)];
        %% We solve the inflated NLEVP
        ev = evs(j);
        evec = [evec; zeros(size(Lambda,1))];
        [ev, evec, resid, normF] = iRefinement(Ftil, ev, evec, u, thres, ...
            delta, ...
            NewtRefs, normF, verbose);
        Y = evec(1:n);
        V = evec(n+1:end);
        
        % The newer eigenpair is surely different from the ones in
        % the current minimal invariant pair, but it could be equal
        % to another one in another minimal invariant pair
        [NewEigenpair, epairIndexNew] = iCheckNewEig(ev, Y, ...
            pairs, tol1);
        pairs.lambdas{epairIndex} = [Lambda V; zeros(1, length(V)) ev];
        pairs.Xs{epairIndex} = [X Y];
        
        % Now we check if the minimality index has changed
        % (=increased). The paper of Effemberger tells us that it
        % increases at most by one.
        sameMinIndex = iCheckIndex(Lambda, X, minIndex);
        if ~sameMinIndex
            pairs.minIndeces(epairIndex) = minIndex+1;
        end
        %If we are very unlucky, we need to repeat this deflation
        % process with another minimal invariant pair, so we update the
        % minimal invariant pair we are considering thanks to the
        % variable 'epairIndex'. Otherwise we set Y as the new
        % eigenvector. We do not need to the same thing for the
        % eigenvalue because we are using the same variable 'ev'.
        if ~NewEigenPair
            epairIndex = epairIndexNew;
        else
            evec = Y;
        end
    end
    
    % We update the eigenvalues, eigenvectors. We do not update the
    % minimal invariant pairs (the variable 'pairs') because we have
    % already done that, either when we entered the deflation cycle or
    % earlier from lines 30 to 35.
    evs(j) = ev;
    evecs(:,j) = evec;
    resids(j) = resid;
end



end

function [newEig, epairIndex] = iCheckNewEig(ev, evec, pairs, tol1)
% This function checks whether the eigenpair computed by the Newton
% refinement is a new one. It goes through the minimal invariant
% pairs. It returns  a boolean 'newEig' and the position of the pair
% 'epairIndex', which is =-1 if 'newEig' is true.

epairIndex = -1;
newEig = true;
for j = 1 : pairs.length
    m  = size(pairs.lambdas{j}, 1);
    %Usually m = 1
    for k = 1 : m
        distEv = abs( pairs.lambdas{j}(k,k) - ev);
        % We check if same eigenvalue
        if distEv < tol1 * abs(pairs.lambdas{j}(k,k))
            % If same eigenvalue, we check if it is the same
            % eigenvector. Given that they
            % may differ by a complex number of norm 1, we need to
            % normalize. We
            % divide by the maximum element. We do not need to normalise
            % because the norm is already 1.
             [~, maxInd] = max(abs(evec));
             distEc = norm(evec - ...
                 pairs.Xs{j}(:,k)*evec(maxInd)/ pairs.Xs{j}(maxInd,k));
            if distEc < tol1 *norm( pairs.Xs{j}(:,k) )
                % We have already computed this eigenpair and it was saved
                % in the jth minimal eigenpair.
                newEig = false;
                epairIndex = j;
                return
            end
        end
    end
end

end

function [ev, evec, resid, normF] = iRefinement(F, ev, evec, u, thres, ...
    delta, NewtRefs, normF, verbose)

for iter = 1 : NewtRefs
    Fev = F(ev);
    if verbose >= 2
        fprintf('NewtRefs: %d,   Orig abs(F(ev)): %2.4e\n', NewtRefs, norm(Fev, 'fro'));
    end
    
    %% random direction to approx the derivative
    v = delta*exp(2*pi*1i*rand(1));
    Fpev = (F(ev + v) - Fev)/v; %secant approximation

    
    
    % Newton step
    evect = Fev\(Fpev * evec);
    ev = ev - u'*evec/(u'*evect);
    evec = evect/norm(evect);
    % End of Newton step
    
    %resid
    normF = max(normF, norm( F(ev), 'fro' ));
    %resid =  norm( F(ev)*evec )/normF;
    resid =  norm( F(ev)*evec );
    
    if verbose >=2
        fprintf('%d normF: %2.4e,  abs(ev): %2.4e   Resid: %2.4e\n', iter, normF, abs(ev),  resid)
    end
    if resid < thres
        if verbose >= 1
            fprintf(['Convergence reached in %d iteration,',  ...
                ' with a resid equal to is %7.1e\n\n'], iter, ...
                resid);
        end
        break
    end
    if iter == NewtRefs
        if verbose >= 1
            warning(['Reached maximum number of Newton refinement %d,',  ...
                ' but the resid for eigenvalue %7.1e is %7.1e\n\n'], NewtRefs, ...
                ev, resid);
        end
    end
end

end

function [A,B,U] = iBuildAuxFun(pairs, epairPos, F)

% Build the functions A(\mu), B(\mu) and U(\mu) as described by
% Effenberger. This is far from the best implementation (aong other
% tweaks, Horner rule, but I think is fine for now
Lambda = pairs.lambdas{epairPos};
X = pairs.Xs{epairPos};
l = pairs.minIndeces(epairPos);
A = @(z) iAA(z, Lambda, X, l);
B = @(z) iBB(z, Lambda, X, l);
U = @(z) iUU(F, z, Lambda, X);
end

function res = iAA(mu, Lambda, X, p)
I = eye(size(Lambda));
res = I;
for j = 1:p
    res = I +  mu*Lambda'*res;
end
res = res*X';
end

function res = iBB(mu, Lambda, X, p)
XTX = X'*X;
res = Lambda'*XTX;
for j = 2:p
    q = 0;
    for k = 0:j-1
        q = q + mu^k*Lambda^(j-1-k);
    end
    res = res + (Lambda^j)'*XTX*q;
end
end

function sameIndex = iCheckIndex(Lambda, X, l)
% We check if the minimality index of the minimal pair (Lambda X) is
% equal to l. We build the matrix [X; X\Lambda; ...; X\Lambda^{-1}]
% and check its rank.
[n,m] = size(X);
minmatrix = zeros(n*l, m);
XLambda = X;
for j = 0:l-1
    minmatrix(1+j*n : n+j*n, :) = XLambda;
    XLambda = XLambda*Lambda;
end
rk = rank(minMatrix);
if rk < m
    sameIndex = false;
else
    sameIndex = true;
end


end

function res = iUU(F, mu, Lambda, X)
% This function builds an approximation of
%
%    U(x) = F(x)X(xI - Lambda)^-1
%
% as proposed by Françoise
I = eye(size(Lambda));
res = F(mu)*X/(mu*I - Lambda);
end

function [IDX, C, D] = iChoiceOfKClusters(n, w, gam, rad, points, ...
    nBadEvs, visual)

% [IDX, C, D] = iChoiceOfKClusters(n, w, gam, rad, points, nBadEvs, visual)
% is an auxiliary function of CONTOURSOLVER. It is used when the recursive
% refinement is called. It divides the bad eigenvalues that need to be
% refined in K clustes of centers C and radii D thanks to the IDX variable.
% If the input "visual" is true, then the number of clusters K is chosen by
% the user at runtime, otherwise it is chosen thanks to the "elbow method"
% heuristic.

badEvs = points(:,1) + 1i*points(:,2);
if visual
    figure(1);
    hold off
    plot(gam +rad*w, '-');
    hold on
    plot(badEvs, '.');
    axis equal;
    grid;
    K = input('Choose a number of clusters:  [3] ');
    if ~((isnumeric(K)) && (isequal(size(K), [1 1])))
        warning('You have not inserted an integer. Setting K = 2.')
        K = 2;
    else
        if ~isequal(K, round(K))
            warning('You have not inserted an integer. Setting K = 2.')
            K = 2;
        end
    end
    [IDX, C, ~, D] = kmeans(points, K);
    plot(C(:,1)+1i*C(:,2), 'o');
else    
    K1 = max(floor(nBadEvs/n),1);
    IDX1 = zeros(nBadEvs, 3);
    C1 = cell(1,3);
    SUMD1 = zeros(max(K1+1,3),3);
    D1 = cell(1,3);
    clusters = max(1,K1-1):max(K1+1,3);
    for j = 1:3
        [IDX1(:,j), C1{j}, SUMD, D1{j}] = kmeans(points, clusters(j));
        SUMD1(1:clusters(j),j) = SUMD;
    end
    err = sum(SUMD1);
    linApprox = (err(3)+err(1))/2;
    plot(err, 'r-o')
    if (err(2) < linApprox*0.75)
        IDX = IDX1(:,2);
        C =  C1{2};
        D = D1{2};
    else
        if  (err(2) > linApprox*1.25)
            IDX = IDX1(:,3);
            C =  C1{3};
            D = D1{3};
        else
            IDX = IDX1(:,1);
            C =  C1{1};
            D = D1{1};
        end
    end
end
end

%% Parsing Inputs

% The purpose of the following auxiliary functions is to parse all the
% inputs chosen by the user and set the defaul parameters of the algorithm.

function [gam, rad, str, rot, n, nc, w, z, M, Vexists, V, Vsize, Fp, ...
    Fpexists, thres, maxRefine, fastUp, onlyIn, ...
    NewtRefs, ref2, nF, maxRecLvl, gk, method, lData, verbose, info] = ...
    iCheck(F, gam, rad, varargin)

if nargin  > 4
    info.err = 'Too many inputs';
    error('myfuns:contourSolver:TooManyInputs', ...
        ['At  most four inputs: the function, the center,' ...
        ' the radius and the options as a structure.']);
end

default_ctrl = 0;
switch nargin
    case 4
        opts = varargin{1};
    case 3
        default_ctrl = 1;
        opts = struct();
    case 2
        default_ctrl = 1;
        opts = struct();
        rad = 1;
    case 1
        default_ctrl = 1;
        opts = struct();
        rad = 1;
        gam = 0;
end

if ~isa(F, 'function_handle')
    info.err = 'F is not a function_handle';
    error('myfuns:contourSolver:NotAFunction', ...
        'F is not a function_handle');
    
end

%% Extract the size of F; check if it is a square matrix
[n1,n2] = size(F(1));
if (n1 ~= n2)
    info.err = 'F(z) is not a square matrix.';
    error('myfuns:contourSolver:NotASquareMatrix', ...
        'F(z) is not a square matrix.');
end
n = n1;

if ~isa(opts, 'struct')
    info.err = 'opts is not struct';
    error('myfuns:contourSolver:NotAStruct', ...
        'opts is not a function_handle');
end

%% Checking center and radius
if ~(isnumeric(gam) && isnumeric(rad) && isscalar(gam) && isscalar(rad))
    info.err = 'The radius or the center are not numbers.';
    error('myfuns:contourSolver:NotANumber', ...
        'The radius or the center are not numbers.');
end

if rad <= 0
    info.err = 'The radius is not strictly positive';
    error('myfuns:contourSolver:NotStrictlyPositive', ...
        'The radius must be a strictly positive number');
end

%% Assigning the default values. Return if 'opts' is not given
[nc, str, rot, w, z, M, Vexists, V, Vsize, Fp, Fpexists, thres, ...
    maxRefine, fastUp, onlyIn, NewtRefs, ref2, nF, maxRecLvl, gk, ...
    method, lData, verbose, info] = iDefault_values(n, gam, rad);

if default_ctrl
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extracting the optional values from opts %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checking the verbosity

if isfield(opts, 'verbose')
    verbose = opts.verbose;
    if ~(verbose == 0 || verbose == 1 || verbose == 2)
        info.err = 'verbose is not 0, 1 or 2';
        error('myfuns:contourSolver:NotValid', ...
            'opts.verbose accepts only 0, 1 and 2 as parameters');
    end
else
    verbose = 0;
end

%% Optional field for Gauss-Kronrod
if isfield(opts, 'GK') && (opts.GK == 1)
    gk.points = [-0.991455371120813 -0.949107912342759 -0.864864423359769 ...
        -0.741531185599394 -0.586087235467691 -0.405845151377397 ...
        -0.207784955007898 0 0.207784955007898 0.405845151377397 ...
        0.586087235467691 0.741531185599394 0.864864423359769 ...
        0.949107912342759 0.991455371120813];
    gk.weights = [0.022935322010529 0.063092092629979 0.104790010322250 ...
        0.140653259715525 0.169004726639267 0.190350578064785 ...
        0.204432940075298 0.209482141084728 0.204432940075298 ...
        0.190350578064785 0.169004726639267 0.140653259715525 ...
        0.104790010322250 0.063092092629979 0.022935322010529];
    gk.weightsGauss = [0.129484966168870, 0.279705391489277, 0.381830050505119, ...
        0.417959183673469, 0.381830050505119, 0.279705391489277, ...
        0.129484966168870];
    if isfield(opts, 'str')
        str = opts.str;
        % Checking the stretch
        if ~(isnumeric(str) && isscalar(str))
            info.err = 'The ratio between  the vertical and the horizontal axis is not a number.';
            error('myfuns:contourSolver:NotANumber', ...
                'The ratio between the vertical and the horizontal axis is not a number.');
        end
        if str <= 0
            info.err = 'The stretch is not strictly positive';
            error('myfuns:contourSolver:NotStrictlyPositive', ...
                'The ratio between the vertical and the horizontal axis must be a strictly positive number');
        end
    else
        str = 1;
    end
    if isfield(opts, 'rot')
        rot = opts.rot;
        % Checking the rotation
        if ~(isnumeric(rot) && isscalar(rot))
            info.err = 'The rotation of the ellipse is not a number.';
            error('myfuns:contourSolver:NotANumber', ...
                'The rotation of the ellipse is not a number.');
        end
    else
        rot = 0;
    end
    
    theta = pi*(gk.points + 1);
    w = cos(theta)*cos(rot) - str*sin(theta)*sin(rot)    +     1i*sin(rot)*cos(theta) + 1i*str*cos(rot)*sin(theta);
    z = gam + rad*w;
    nc = 15;
else
    if isfield(opts, 'nc')
        nc = opts.nc;
        if ~(mod(nc,1) == 0 && isscalar(nc))
            info.err = 'nc is not a number';
            error('myfuns:contourSolver:NotANumber', ...
                'nc is not a number');
        end
        if nc <= 0
            info.err = 'The number of integration points is negative';
            error('myfuns:contourSolver:NotStrictlyPositive', ...
                'The number of integration points must be strictly positive.');
        end
    end
    if isfield(opts, 'str')
        str = opts.str;
        % Checking the stretch
        if ~(isnumeric(str) && isscalar(str))
            info.err = 'The ratio between  the vertical and the horizontal axis is not a number.';
            error('myfuns:contourSolver:NotANumber', ...
                'The ratio between the vertical and the horizontal axis is not a number.');
        end
        if str <= 0
            info.err = 'The stretch is not strictly positive';
            error('myfuns:contourSolver:NotStrictlyPositive', ...
                'The ratio between the vertical and the horizontal axis must be a strictly positive number');
        end
    else
        str = 1;
    end
    if isfield(opts, 'rot')
        rot = opts.rot;
        % Checking the rotation
        if ~(isnumeric(rot) && isscalar(rot))
            info.err = 'The rotation of the ellipse is not a number.';
            error('myfuns:contourSolver:NotANumber', ...
                'The rotation of the ellipse is not a number.');
        end
    else
        rot = 0;
    end
    theta = 2*pi*(1:nc)/nc;
    w = cos(theta)*cos(rot) - str*sin(theta)*sin(rot)    +     1i*sin(rot)*cos(theta) + 1i*str*cos(rot)*sin(theta);
    z = gam + rad*w;
end

if isfield(opts, 'M')
    M = opts.M;
    if ~(mod(M,1) == 0 && isscalar(M))
        info.err = 'M is not a number';
        error('myfuns:contourSolver:NotANumber', ...
            'M is not a number.');
    end
    if (M <= 0)
        info.err = 'The number of moments is negative.';
        error('myfuns:contourSolver:NotStrictlyPositive', ...
            'The number of moments must be a strictly positive.');
    end
end

% If the number of moments M is greater than the number of quadrature points
% for the trapezoidal rule nc, than A_{M+1} would be equal to A_{1} due to
% its approximation. The following check avoids this issue.

if ~(isfield(opts, 'GK') && (opts.GK == 1)) && M > nc
    info.err = ['The number of moments is greater than the number '...
        'of trapezoidal quadrature points'];
    error(['myfuns:contourSolver: ', ...
        'The number of moments is greater than the number ', ...
        'of trapezoidal quadrature points']);
end

%% Standard Hankel or Loewner interpretation
if isfield(opts, 'method')
    method = opts.method;
    if ~(isequal(method, 'std') || isequal(method, 'loewner'))
        info.err = 'opts.method only accepts "std" or "loewner".';
        error('myfuns:contourSolver:Method', ...
            'opts.method only accepts "std" and "loewner" as parameters.');
    end
end

%% Extraxt the starting projection space
if isfield(opts, 'V')
    V = opts.V;
    if ~(isnumeric(V))
        info.err = 'V is not a matrix.';
        error('myfuns:contourSolver:NotNumeric', ...
            'The projection space V is not a matrix.');
    end
    [Vcols, Vsize] = size(V);
    if Vcols ~= n
        info.err = 'Mismatching dimensions.';
        error('myfuns:contourSolver:DimensionMismatch', ...
            'Mismatching dimension between F(z) and V.');
    end
    Vexists = true;
else
    Vexists = false;
end

%% Extracting the optional size of the projection space

if ~Vexists
    if isfield(opts, 'Vsize')
        if ~(isscalar(opts.Vsize))
            info.err = 'Vsize is not a number';
            error('myfuns:contourSolver:NotANumber', ...
                'Vsize is not a number.');
        else
            if (mod(opts.Vsize, 1) ~= 0 || opts.Vsize <= 0)
                info.err = ['The dimension of the projection space is not ' ...
                    'a strictly positive integer.'];
                error('myfuns:contourSolver:NotStrictlyPositive', ...
                    ['The dimension of the subspace must be ' ...
                    'a strictly positive integer.']);
            end
            if opts.Vsize > n && isequal(method, 'std')
                warning(['Vsize must be <= the size of the F(z).' ...
                    ' Automatically set Vsize = %d'], n);
            end
            Vsize = min(n, round(opts.Vsize));
            V = iBuildV(n, Vsize);
        end
    end
    
    
    %% V is given
else
    if (isfield(opts, 'Vsize'))
        if verbose
            warning(['Ignoring the parameter opts.Vsize because the user',  ...
                'already gave V. Using instead Vsize = size(V,2).']);
        end
    end
end

%% Extracting the Loewner points and directions, if they exist
if isfield(opts, 'lData')
    lData = opts.lData;
    lData.exists = 1;
    %Check the number of vectors is the same number of points
    if length(lData.points) ~= size(lData.dir, 2)
        info.err = ['The number of interpolation points is different from',...
            ' the number of directions.'];
        error('myfuns:contourSolver:wrongInterpolationData', ...
            ['The number of interpolation points is different from' ...
            ' the number of directions.']);
    end
else
    lData.exists = 0;
end

%% Eigenvalues in or out?

if isfield(opts, 'onlyIn')
    onlyIn = opts.onlyIn;
    if ~(onlyIn == 1 || onlyIn == 0)
        error('myfuns:contourSolver:NotABool', ...
            'onlyIn is not a boolean.');
    end
end

%% Derivative of F

if isfield(opts, 'Fp') && ~isempty(opts.Fp)
    Fp = opts.Fp;
    if ~isa(Fp, 'function_handle')
        info.err = 'Fp is not a function_handle';
        error('myfuns:contourSolver:NotAFunction', ...
            'Fp is not a function_handle');
    end
    Fpexists = true;
else
    Fpexists = false;
    Fp = [];
end

%% Threshold for the eigenvalues

if isfield(opts, 'thres')
    thres = opts.thres;
    if ~isnumeric(thres)
        info.err = 'thres is not a numbe.r';
        error('myfuns:contourSolver:NotANumber', ...
            'thres is not a number.');
        
    end
end

%% How many probing refinements

if isfield(opts, 'maxRefine')
    maxRefine = opts.maxRefine;
    if ~isnumeric(maxRefine)
        info.err = 'maxRefine is not a number';
        error('myfuns:contourSolver:maxRefine', ...
            'maxRefine is not a number.');
    end
end

%% Fast update

% This will work in a future release
if isfield(opts, 'fastUp')
    fastUp = opts.fastUp;
    if fastUp == 1
        %% check if colqr actually exists, otherwise throw an error
        if (exist('colqr') ~= 3)
            info.err = 'fastUp = 1 without colqr';
            error('myfuns:contourSolver:fastUp', ...
                ['You can not set fastUp = 1 if you do not have ' ...
                'colqr.']);
            
        end
    elseif fastUp == 0
        if verbose > 0
            fprintf(['Not using a fast rank update to find the optimal', ...
                ' size of the matrix V.\n']);
        end
    else
        info.err = 'fastUp is not 1 neither 0';
        error('myfuns:contourSolver:fastUp', ...
            'fastUp is not 1 neither 0.');
    end
end

%% How many Newton refinements

if isfield(opts, 'NewtRefs')
    NewtRefs = opts.NewtRefs;
    if ~isnumeric(NewtRefs)
        info.err = 'NewtRefs is not a number';
        error('myfuns:contourSolver:NewtRefs', ...
            'NewtRefs is not a number.');
    end
end

%% Which refinement

if isfield(opts, 'ref2')
    ref2 = opts.ref2;
    if ~(ref2 == 0 || ref2 == 1 || ref2 == 2 || ref2 == 3 || ref2 == 4)
        info.err = 'ref2 is not an accepted value';
        error('myfuns:contourSolver:ref2', ...
            'ref2 is not 1 neither 0.');
    end
end

%% Recursion level

if isfield(opts, 'maxRecLvl')
    maxRecLvl = opts.maxRecLvl;
    if ~isnumeric(maxRecLvl)
        info.err = 'maxRecLvl is not a number';
        error('myfuns:contourSolver:MaxRecLvl', ...
            'MaxRecLvl is not a number.');
    end
end

%% ||F|| norm

if isfield(opts, 'nF')
    if isnumeric(opts.nF) && opts.nF >= 0
        nF = opts.nF;
    else
        info.err = 'opts.nF is not a non-negative real number.';
        error('myfuns:contourSolver:nF', ...
            'opts.nF must be a non-negative real number.');
    end
    
end

%% Setting the info
info.gam = gam;
info.rad = rad;
info.nc = nc;
info.str = str;
info.rot = rot;
info.w = w;
info.z = z;
info.M = M;
info.Vexists = Vexists;
info.V = V;
info.Vsize = Vsize;
info.Fpexists = Fpexists;
info.Fp = Fp;
info.thres = thres;
info.maxRefine = maxRefine;
info.fastUp = fastUp;
info.nc = nc;
info.method = method;
info.verbose = verbose;
info.onlyIn = onlyIn;
info.NewtRefs = NewtRefs;

end


function [nc, str, rot, w, z, M, Vexists, V, Vsize, Fp, Fpexists, ...
    thres, maxRefine, fastUp,  onlyIn, NewtRefs, ref2, nF, maxRecLvl, ...
    gk, method, lData, verbose, info] = iDefault_values(n, gam, rad)

% [nc, w, z, M, Vexists, V, Vsize, Fp, Fpexists, ...
% thres, maxRefine, fastUp,  onlyIn, verbose, info] = ...
% DEFAULT_VALUES(F,n)
% is an auxiliary function of  CONTOURSOLVER2. It is used to
% set the default values when not given by the user.
%
% See also CONTOURSOLVER.

str = 1;
rot = 0;
nc = 2*max(round(16*2^log10(rad)), 2);
theta = 2*pi*(1:nc)/nc;
w = exp(1i*theta);
z = gam + rad*w;
M = 2;
Vexists = false;
Vsize = 1;
V = iBuildV(n,Vsize);
Fp = [];
Fpexists = false;
thres = 1e-11;
maxRefine = 1;
fastUp = 0;
onlyIn = true;
NewtRefs = 5;
ref2 = 2;
nF = 0;
maxRecLvl = 1;
gk = struct();
method = 'std';
verbose = 0;

%% Default values for Loewner method
lData.dir = iBuildV(n,2*Vsize);
lData.dir = lData.dir./vecnorm(lData.dir);
lData.points = gam + rad*eps^(-1/nc)*w(1:2*Vsize);
%% Default info
info.gam = gam;
info.rad = rad;
info.nc = nc;
info.w = w;
info.z = z;
info.M = M;
info.Vexists = Vexists;
info.V = V;
info.Vsize = Vsize;
info.Fpexists = Fpexists;
info.Fp = Fp;
info.thres = thres;
info.maxRefine = maxRefine;
info.fastUp = fastUp;
info.onlyIn = onlyIn;
info.NewtRefs = NewtRefs;
info.verbose = verbose;
info.ref2 = ref2;
info.normF = nF;
info.method = method;
info.err = '';

end

