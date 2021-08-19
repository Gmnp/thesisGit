function fig = plot_newton_polygon(fi, fa, I, V, a, b, varargin)
%PLOT_NEWTON_POLYGON 

% if ~exist('nnodes', 'var')
%     nnodes = 10;
% end

nnodesl = 10;
nnodesr = 10;

if nargin == 7
    nnodesl = varargin{1};
    nnodesr = nnodesl;
else % nargin = 8
    nnodesl = varargin{1};
    nnodesr = varargin{2};
end

Ir = []; Vr = [];
for j = 0 : nnodesr - 1
    Ir = [ Ir, fi(b+j) ];
    Vr = [ Vr, fa(b+j) ];
end

Il = []; Vl = [];
for j = 0 : nnodesl - 1
    Il = [ fi(a-j), Il ];
    Vl = [ fa(a-j), Vl ];
end

I = [ Il, I, Ir ];
V = [ Vl, V, Vr ];

fig = plot(I, V, '*-');

xlabel('$j$', 'Interpreter', 'latex'); ylabel('$\log b_j$', 'Interpreter', 'latex');


end

