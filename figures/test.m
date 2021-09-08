mydefaults % call this first to initialize font size and line thickness

% create some plot
plot(randn(50,3))

mypdf('filename',0.7,2) % save images

% 1st parameter: file name without extension
% 2nd parameter: width to height ratio (1 would be almost square)
% 3rd parameter: scaling, increase if font size is too big