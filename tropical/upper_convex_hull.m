function [I, V] = upper_convex_hull(I, V)
%UPPER_CONVEX_HULL2 

x = min(V) - 1;
k = convhull([ min(I), I, max(I)  ], [ x, V, x ]);
k = sort(k);
k = k(3:end-1) - 1;
I = I(k); V = V(k);

end

