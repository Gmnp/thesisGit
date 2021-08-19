function [I, V, a, b] = newton_polygon(fi, fa, ip, p)
%NEWTON_POLYGON Update an infinite Newton Polygon with a polynomial. 
%
% This function takes as input two functions that, for a given integer I,
% return the index fi(I) and the value fa(I) of a vertex of the Newton
% polygon at that index, and a vector P. It is expected that fi(I) has the
% same sign of I. 
%
% The returned vectors I, V are a set of indices such that the Newton
% Polygon of f(x) + p(x) is described by the vertices
%
%   I(1) ... I(k)
%  
% with value V(1), ..., V(k), and the infinite number of points of
% coordinates (fi(I), fa(I)) for the indices I <= A or I >= B. It can
% happen that A or B are equal to -inf and inf, respectively. 

% Step 1: we find the extrema a and b such that the indices ip are in the
%         interval fi(A) < ... < fi(B); this way, we can update the indices
%         of the vector P including the ones in FA(I) for these indices. 

imax = max(ip); imin = min(ip);

b = 0; fj = fi(b);
while fj <= imax
    b = b + 1;
    fj = fi(b);
end

a = 0; fj = fi(a);
while fj >= imin
    a = a - 1;
    fj = fi(a);
end

% Step 2: Merge the two vectors of points (ip, p) and (fi, fa) in the 
%         interval a .. b. We do this by first computing the indices of the
%         original tropcal function in a .. b, and then we scan the arrays
%         entry by entry and take the maximum in common positions. 
ni = []; na = [];
for i = a+1 : b-1
    ni = [ni, fi(i)];
    na = [na, fa(i)];
end

jt = 1; jp = 1;
I = []; V = [];
while jp <= length(p) && jt <= length(ni)
    if ip(jp) < ni(jt)
        I = [ I, ip(jp) ];
        V = [ V, p(jp) ];
        jp = jp + 1;
    elseif ni(jt) < ip(jp)
        I = [ I, ni(jt) ];
        V = [ V, na(jt) ];
        jt = jt + 1;
    else
        % In this case the indices are equal, we take the maximum
        I = [ I, ni(jt) ];
        V = [ V, max(na(jt), p(jp)) ];
        jt = jt + 1; jp = jp + 1;
    end
end

% Step 3: we now compute the upper convex hull of the points in [I, V], and
% we drop the indices that are not part of it. 
[I, V] = upper_convex_hull(I, V);

% Step 4: Merge the convex hull defined by I, V with the rest of the
% points; we start by the ones on the right. We do this as follows:
%
%  -) We select as first candidate the semgent with indices I(end) -- fi(b)
%  -) We compute its slope
%  -) We check if it is convex with the pieces on the left and right; if
%     not, we increment / decrement the indices, and continue until
%     convergence
%
% FIXME: We need to handle the case with only one point. 
hull_converged = false;

while ~hull_converged
    hull_converged = true;
    
    slope = ( fa(b) - V(end) ) / ( fi(b) - I(end) );
    
    if length(V) > 1
        lslope = ( V(end) - V(end-1) ) / ( I(end) - I(end-1) );
    else
        lslope = inf;
    end
    
    rslope = (fa(b+1) - fa(b)) / (fi(b+1) - fi(b));
    
    if lslope < slope
        I = I(1:end-1);
        V = V(1:end-1);
        hull_converged = false;
    end
    
    if rslope > slope
        b = b + 1;
        hull_converged = false;
    end
end


% Step 4.5: The same as Step 4, but on the left
hull_converged = false;

while ~hull_converged
    hull_converged = true;
    
    slope = ( V(1) - fa(a) ) / ( I(1) - fi(a) );
    
    if length(V) > 1
        rslope = ( V(2) - V(1) ) / ( I(2) - I(1) );
    else
        rslope = -inf;
    end
    
    lslope = (fa(a) - fa(a-1)) / (fi(a) - fi(a-1));
    
    if rslope > slope
        I = I(2:end);
        V = V(2:end);
        hull_converged = false;
    end
    
    if lslope < slope
        a = a - 1;
        hull_converged = false;
    end
end

end