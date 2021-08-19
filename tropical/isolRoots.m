nPoints = 100;
x = linspace(1, 14, nPoints);
y = linspace(0.001, 1, nPoints);
z = [ linspace(0.001, 0.17, nPoints), linspace(0.17, 14, nPoints)];
z = linspace(0.001, 14, 4*nPoints);
% To get the right colors
aux = lines;

figure(1)
hold off
plot([y, x], 3*ones(1,2*nPoints), 'Color', aux(1,:))
hold on
plot([y x], 0.5*[y x].^(1), 'Color', aux(2,:))
plot([y x], 0.2*[y x].^(-1), '-.', 'Color', aux(2,:))
plot([y x], [y x].^(2)/16, 'Color', aux(3,:))
plot([y x], [y x].^(-2)/2, '-.', 'Color', aux(3,:))
g2 = max([3*ones(length(z)); 0.5*z.^(1); 0.2*z.^(-1); z.^(2)/16; z.^(-2)/2]); % we take the max
plot(z, g2, '--k', 'LineWidth', 2.5)
legend({'$b_0$', '$b_1x$', '$b_{-1}x^{-1}$', '$b_2x^2$', '$b_{-2}x^{-2}$', '$g_2(x)$'}, 'Interpreter', 'latex', 'Location', 'best', 'NumColumns',2)
axis([0, 10, 0,6])


 figure(2)
 hold off
plot([y, x], 3*ones(1,2*nPoints), 'Color', aux(1,:))
hold on
plot([y x], 0.5*[y x].^(1), 'Color', aux(2,:))
plot([y x], 0.2*[y x].^(-1), '-.', 'Color', aux(2,:))
plot([y x], [y x].^(2)/12, 'Color', aux(3,:))
plot([y x], [y x].^(-2)/2, '-.', 'Color', aux(3,:))
g2 = max([3*ones(length(z)); 0.5*z.^(1); 0.2*z.^(-1); z.^(2)/12; z.^(-2)/2]); % we take the max
plot(z, g2, '--k', 'LineWidth', 2.5)
legend({'$b_0$', '$b_1x$', '$b_{-1}x^{-1}$', '$b_2x^2$', '$b_{-2}x^{-2}$', '$g_2(x)$'}, 'Interpreter', 'latex', 'Location', 'best', 'NumColumns',2)
% hold off
% plot([y x], 3*ones(2*nPoints), 'Color', aux(1,:))
% hold on
% plot(x, 0.5*x.^(1), 'Color', aux(2,:))
% plot(y, 0.5*y.^(-1), '--', 'Color', aux(2,:))
% plot(x, x.^(2)/12, 'Color', aux(3,:))
% plot(y, y.^(-2)/12, '--', 'Color', aux(3,:))
% g2 = max([3*ones(length(z)); 0.5*z.^(1); 0.5*z.^(-1); z.^(2)/12; z.^(-2)/12]); % we take the max
% plot(z, g2, '--k')
axis([0, 10, 0,6])