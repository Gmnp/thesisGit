function [evs, evsClean] = contourRIM(F, gam, xhalfedge, yhalfedge, S3, verbose, varargin)

if nargin < 4
    error('This function needs at least 4 inputs')
end
if nargin < 5
S3 = 3;
end
if nargin < 6
verbose = 0;   
end
if nargin > 6
    firstcall = varargin{1};
else
    firstcall = 1;
end

evsClean = [];
Fp = [];
Fpexists = 0;
n = size(F(1),1);
lx = xhalfedge;
ly = yhalfedge;

sqSize = sqrt(lx^2+ly^2);





% Check integration
nc = 2;
% Points of the rectangle
z = linspace(-lx-1i*ly, lx-1i*ly, nc+1);
z = horzcat(z(1:end-1), linspace(lx-1i*ly, lx+1i*ly, nc+1));
z = horzcat(z(1:end-1), linspace(lx+1i*ly, -lx+1i*ly, nc+1));
z = horzcat(z(1:end-1), linspace(-lx+1i*ly, -lx-1i*ly, nc+1));
z = gam + z(1:end-1);

% Points of the normalized rectangle
w = (z - gam)/lx;

rEg = ly/lx; % rationEdges

%weights
weights = [1-1i*rEg, 2*ones(1,nc-1), 1+1i*rEg, 2i*rEg*ones(1,nc-1)];
weights = [weights -weights]; %*lx/(2i*pi*nc);
% We multiply by lx/(2i*pi*nc) at the end of the cycle of contourCountRIM
% to diminish num errors


lupoints = cell(1, 4*nc);
for k = 1:4*nc
    lupoints{k} = decomposition(F(z(k)));
end

  hold on
  plot([z z(1)], 'Color', lines(1), 'Linewidth', 1.5)

% Indicators
% Indicator based on counting eigenvalues
nevs = contourCountRIM(F, n, 4*nc, weights, lx,  z, lupoints, Fp, Fpexists, verbose);
% original indicator by Huang
%ind = iRIMindicator(n, 4*nc, weights, lx, lupoints);
clear lupoints
threshInd = 0.1;
thresh = 5*1e-2;
if abs(nevs) > threshInd

    if sqSize < thresh
        evs = gam;
    else
        evs = [];
        if S3 == 2
            lxn = lx/2;
            lyn = ly/2;
            gams = gam+[lx/2+1i*ly/2, -lx/2+1i*ly/2, -lx/2-1i*ly/2, lx/2-1i*ly/2];
        else
            lxn = lx/3;
            lyn = ly/3;
            gams = gam + [-2*lxn+2i*lyn, 2i*lyn, 2*lxn+2i*lyn, -2*lxn, 0, 2*lxn, -2*lxn-2i*lyn, -2i*lyn, 2*lxn-2i*lyn];
        end
        for j = 1:length(gams)
            evsRec = contourRIM(F, gams(j), lxn, lyn, S3, verbose, 0);
            evs = [evs; evsRec];
        end
    end
else
    evs = [];
end

if firstcall
   evsClean=evs; 
   pairs = abs(evs(1:end-1)-evs(2:end)) <= 2*thresh;
   indexPairs = find(pairs);
   evsClean(indexPairs) = (evs(indexPairs)+ evs(indexPairs+1))/2;
   evsClean(indexPairs+1) = [];  
end

end




function ind = iRIMindicator(n, nc, weights, lx, lupoints)

rng(42);
p = randn(n,1);
newp = zeros(n,1);
newp2 = newp;
% newp = Ind(p)
for k = 1:nc
newp = newp + weights(k)*lupoints{k}\p;
end
newp = 4*newp*lx/(2i*pi*nc);
Alpha=1;
newp = Alpha*newp;
% ind = || Ind(Ind(p)) ||_2
for k = 1:nc
  newp2 = newp2  + weights(k)*lupoints{k}\newp;
end
ind = 4*Alpha*norm(newp2)*lx/(2*pi*nc) / norm(newp);
end