%% Function f(z) = exp(z) + exp(1/z) - 1, with some coefficients changed. 

close all;


fi = @(i) i;
fa = @(i) sum(arrayfun(@(j) -log(j), 1 : abs(i)));

iq  = [ -9 , -3 , 1, 2, 3, 4, 5 ];
q   = exp([ 6, 12, 0,  2, -10, -14, -20 ]);

[II, VV, a, b] = newton_polygon(fi, fa, iq, log(q + 1./factorial(abs(iq))));

% To get the right colors
colors = lines;

%fig1 = plot_newton_polygon(fi, fa, [], [], 0, 0, 50);
fig1 = plot_newton_polygon(fi, fa, [], [], 0, 0, 42);
fig1.LineWidth = 1.2;
hold on;
%plot_newton_polygon(fi, fa, II, VV, a, b);
fig11 = plot_newton_polygon(fi, fa, II, VV, a, b, 4, 38);
fig11.LineWidth = 1.2;
fig12 = plot(iq, log(q + 1./factorial(abs(iq))), 'o', 'Color', colors(2,:), 'LineWidth', 1.2, 'MarkerSize',9);

% We now find the slopes and their multiplicity
% Construct the nodes of the Newton polygon
iA = [ arrayfun(@(j) fi(j), a-100 : a-1), II, arrayfun(@(j) fi(j), b+1 : b+10) ];
A  = [ arrayfun(@(j) fa(j), a-100 : a-1), VV, arrayfun(@(j) fa(j), b+1 : b+10) ];

% And the slopes and multiplicities
A = exp(- (A(2:end) - A(1:end-1)) ./ (iA(2:end) - iA(1:end-1)) );
iA = iA(2:end) - iA(1:end-1);

[ir,er,nr] = annuli(iA, A);

% find the roots?
ff = @(x) exp(x) + exp(1/x) - 1 + sum(q .* (x.^(iq)));
fd = @(x) exp(x) - exp(1/x) ./ (x.^2) + sum((iq .* q) .* (x.^(iq - 1)));
r = [];
t = linspace(-5, 5, 101);
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

t = linspace(-ir(1), ir(1), 201);
for i = 1 : length(t)
    for j = 1 : length(t)
        x = t(j) + 1i * t(i);
        newt_corr = inf;
        k = 1;
        while abs(newt_corr) > 1e-15 && k < 1000
            newt_corr = ff(x) / fd(x);
            x = x - newt_corr;
            k = k + 1;
        end
        
        if (isempty(r) || min(abs(x - r)) > 1e-6) && abs(x) < 10 && k < 1000
            r = [r, x];
        end
    end
end

% find the poles?
poles = [];
% t = linspace(-ir(1), ir(1), 101);
% for i = 1 : length(t)
%     for j = 1 : length(t)
%         x = t(j) + 1i * t(i);
%         newt_corr = inf;
%         k = 1;
%         while abs(newt_corr) > 1e-15 && k < 1000
%             newt_corr = -fd(x) / ff(x);
%             x = x - newt_corr;
%             k = k + 1;
%         end
%         
%         if (isempty(poles) || min(abs(x - poles)) > 1e-6) && abs(x) < 10 && k < 1000
%             poles = [poles, x];
%         end
%     end
% end


fig2 = figure;
% we only plot the roots outside the first disk
r1 = r(abs(r)> ir(1));
plot(r1,  '.', 'Color', colors(2,:), 'MarkerSize',12); hold on;
hold on;
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
axis equal;
axis([ -er(end)*1.05, +er(end)*1.05, -er(end)*1.05, +er(end)*1.05]);
xlabel('$\Re(z)$', 'Interpreter', 'latex');
ylabel('$\Im(z)$', 'Interpreter', 'latex');


% Save the images
saveas(fig1, 'example2New_1.png');
saveas(fig2, 'example2New_2.png');

set(0, 'CurrentFigure', fig2);
%axis([ -0.8, 0.8, -0.8, 0.8 ]);
%saveas(fig2, 'example2New_zoom.png');