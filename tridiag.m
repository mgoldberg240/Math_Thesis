% Created by Matt Goldberg on Jan 20th 2018
% solve a tridiagonal linear system Ax = f for vector x
% A = tridiag(l,d,u)
% Algorithm 2.6 from Epperson Numerical Analysis textbook (p. 80)
function f = tridiag(l,d,u,f)
!
! Elimination Stage
!
n = length(d);
for i = 2:n
    d(i) = d(i) - u(i-1)*l(i)/d(i-1);
    f(i) = f(i) - f(i-1)*l(i)/f(i-1);
end

!
! Backsolve Stage
!

x(n) = f(n)/d(n); % bottom row is a special case
for i = n-1:-1:1
    x(i) = (f(i) -u(i)*x(i+1))/d(i);
end

