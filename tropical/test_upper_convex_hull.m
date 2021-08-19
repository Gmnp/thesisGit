figure;

imax = 30;
n = 10;

I = sort(unique(randi(imax, 1, n)));
V = rand(1, length(I));

[II, VV] = upper_convex_hull(I, V);

plot(I, V, 'b-');
hold on;
plot(II, VV, 'r*-');