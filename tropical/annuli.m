function [ir, er, nr] = annuli(ip, p)
%ANNULI Construct exclusion annuli with internal and external radii, and 
% the number of roots minus poles inside. 

ir = []; 
er = []; 
nr = [];

% Compute delta
delta = p(1:end-1)  ./ p(2:end);

for j = 1 : length(delta)
    d = delta(j);
    if d < 1/9
        % Compute f, g
        rr = roots([ 1, -(2+(1-d)./(2*d)), 1./d ]);
        f = min(rr);
        g = max(rr);
        
        ir = [ ir, f * p(j) ];
        er = [ er, g * p(j) ];
        nr = [ nr, ip(j) ];
    end
end

end

