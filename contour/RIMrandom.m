rng(42)
A = rand(20,20);
I = eye(20);
F = @(z) z*I-A;

% Remember to uncomment the plotting if measuring the time
aux = lines(2); % this allows us to have the usual second color of matlab
aux = aux(2,:);
figure(1)
hold off
tic;[evs, evsClean] = contourRIM(F, 0, 1.5,1.5,3,0);toc
plot(evs, '.', 'Color', aux, 'MarkerSize', 15)
mypdf('RIM9', 0.7,2)
figure(2)
hold off
tic;[evs1, evsClean1] = contourRIM(F, 0, 1.5,1.5,2,0);toc
plot(evs1, '.', 'Color', aux, 'MarkerSize', 15)
mypdf('RIM4', 0.7,2)

fprintf('S = 4: #Evs %d,  #(cleaned evs) %d\n', length(evs1), length(evsClean1))
fprintf('S = 9: #Evs %d,  #(cleaned evs) %d\n', length(evs), length(evsClean))