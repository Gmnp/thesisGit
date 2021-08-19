[~,~, F] = nlevp('butterfly');

% Remember to comment the plotting if measuring the time
aux = lines(2); % this allows us to have the usual second color of matlab
aux = aux(2,:);
figure(1)
hold off
tic;[evs, evsClean] = contourRIM(F, 0, 2,2,3,0);toc
butt1 = plot(evs, '.', 'Color', aux, 'MarkerSize', 10)
mypdf('butt9', 0.7,2)
figure(2)
hold off
tic;[evs1, evsClean1] = contourRIM(F, 0,2, 2, 2,0);toc
butt2 = plot(evs1, '.', 'Color', aux, 'MarkerSize', 10)
mypdf('butt4', 0.7,2)

fprintf('S = 4: #Evs %d,  #(cleaned evs) %d\n', length(evs1), length(evsClean1))
fprintf('S = 9: #Evs %d,  #(cleaned evs) %d\n', length(evs), length(evsClean))