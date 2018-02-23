% Created by Matt Goldberg on Jan 20th 2018
% solve a tridiagonal linear system Ax = f for vector x
% A = tridiag(l,d,u); where l = [0 l2 ... ln], u = [u1 u2 ... u_n-1 0]
% Algorithm 2.6 from Epperson Numerical Analysis textbook (p. 80)
function x = tridiag_solve(l,d,u,f)
% --- Elimination stage
n = length(f);

% del(1) = d(1)
% LHS d below is delta in text
% LHS f below is g in text
for i = 2:n
    d(i) = d(i) - u(i-1)*l(i)/d(i-1);
    f(i) = f(i) - f(i-1)*l(i)/d(i-1);
end

% Backsolve stage

x(n) = f(n)/d(n); % bottom row is a special case
for i = n-1:-1:1
    x(i) = (f(i) - u(i)*x(i+1))/d(i);
end

end