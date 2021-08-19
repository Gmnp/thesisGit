close all;

% We define the functions for the exponential; MATLAB makes it particularly
% horrible to write functione defined piecewise inline ...
fi_tmp = { ...
    @(i) -inf, ... % i <= -1
    @(i) i, ... % i >= 0
};

fi = @(i) fi_tmp{ (i <= -1) * (1) + (i >= 0) * 2 }(i);

fa_tmp = { ...
    @(i) -inf, ... % i <= 1 
    @(i) 0, ... % i == 0
    @(i) sum(arrayfun(@(j) -log(j), 1 : i)) ... % i >= 1;
};

fa = @(i) fa_tmp{ ...
    (i <= -1) * 1 + ...
    (i == 0)  * 2 + ...
    (i >= 1)  * 3}(i);

ip = [ 0, 1, 2, 3, 4, 5, 6 ];
% p = [ 3, 4, -inf, 8, -2, 0, 11 ];
%p = [ 0, 2.5, -inf, 2.5, -10, -4, -14 ];
%pp = exp(p(end:-1:1)) - 1 ./ factorial(length(p)-1:-1:0);

%pp = [ 6e-4, 2.2, 1.5e-3, 73, 0, 12, 1 ];
pp = [ -2e-3, 1e-3, -4.e-02, 12, -1/5, 12, 0 ];
p = log(abs(pp(end:-1:1) + 1 ./ factorial(0:length(pp)-1)));

[II, VV, a, b] = newton_polygon(fi ,fa, ip, p);

% To get the right colors
colors = lines;

%fig1 = plot_newton_polygon(fi, fa, [], [], 0, 0, 15);
fig1 = plot_newton_polygon(fi, fa, [], [], 0, 0, 15);
fig1.LineWidth = 1.2;
hold on;
%plot_newton_polygon(fi, fa, II, VV, a, b);
fig11 = plot_newton_polygon(fi, fa, II, VV, a, b,4);
fig1.LineWidth = 1.2;
plot(ip, p, 'o', 'Color', colors(2,:), 'LineWidth', 1.2, 'MarkerSize',9);
% Construct the nodes of the Newton polygon
ia = [ II, arrayfun(@(j) fi(j), b+1 : b+10) ];
a  = [ VV, arrayfun(@(j) fa(j), b+1 : b+10) ];

% And the slopes and multiplicities
a = exp(- (a(2:end) - a(1:end-1)) ./ (ia(2:end) - ia(1:end-1)) );
ia = ia(2:end) - ia(1:end-1);

[ir,er,nr] = annuli(ia, a);

% We draw the special case of the exclusion disc at zero.
ir = [0, ir];
er = [a(1) / 2, er]; % 2 = 1 + c

% Find some roots, computing them on a grid. 
% t = linspace(-30, 30, 101); [xx, yy] = meshgrid(t, 1i * t);
% r = [];
% for i = 1 : length(t)
%     for j = 1 : length(t)
%         x = fsolve(@(z) exp(z) + polyval(pp, z), xx(i) + yy(j));
%         if isempty(r) || min(abs(x - r)) > 1e-3
%             r = [r, x];
%         end
%     end
% end

% Chebfun works quite well for finding some rough estimates, that we then 
% refine with three Newton's iterations. 
% f = chebfun2(@(x,y) exp(x + 1i * y) + polyval(pp, x + 1i * y), [-30 30 -30 30]);
% r = roots(real(f), imag(f), 'ms');
ff = @(x) ( exp(x) + polyval(pp, x) );
fd = @(x) ( exp(x) + polyval(polyder(pp),x) );

r = [];
t = linspace(-6, 6, 101);
for i = 1 : length(t)
    for j = 1 : length(t)
        x = t(j) + 1i * t(i);
        newt_corr = inf;
        
        while abs(newt_corr) > 1e-8
            newt_corr = ff(x) / fd(x);
            x = x - newt_corr;
        end
        
        if (isempty(r) || min(abs(x - r)) > 1e-1) && abs(x) < 10
            r = [r, x];
        end
    end
end


fig2 = figure;
plot(r,  '.', 'Color', colors(2,:), 'MarkerSize', 11); hold on;
t = linspace(0, 2 * pi, 100);
for i = 1 : length(ir)
    x = ir(i) * cos(t);
    y = ir(i) * sin(t);
    X = er(i) * cos(t);
    Y = er(i) * sin(t);
    plot(x, y, '-', 'Color', colors(1,:), 'LineWidth', 1);
    plot(X, Y, '-', 'Color', colors(1,:), 'LineWidth', 1);
    patch([x X],[y Y], colors(1,:),'linestyle','non','facealph',.3);
end
axis equal
axis([ -er(end)*1.05, +er(end)*1.05, -er(end)*1.05, +er(end)*1.05]);
xlabel('$\Re(z)$', 'Interpreter', 'latex');
ylabel('$\Im(z)$', 'Interpreter', 'latex');

% Save the images
saveas(fig1, 'example1_1.png');
saveas(fig2, 'example1_2.png');
