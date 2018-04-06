
nx = 32;
nt = 32; % N in Epperson

% --- x grid
x0 = 0;
x = linspace(x0,xmax,nx)';
h = xmax/(nx-1);
% --- t grid
t0 = 0;
t = linspace(t0,tmax,nt);
dt = tmax/(nt-1);

cnodes = cos((2*x-1)*pi/(2*xmax)); % chebyshev nodes
T = @(n,cnodes) cos(n.*acos(cnodes));

C_even = T(2:2:nt,cnodes)-T(0,cnodes); % for even n
C_odd = T(1:2:nt+1,cnodes)-T(1,cnodes); % for odd n
C=[C_even; C_odd]; C = C(:)' % combine even and odd C arrays