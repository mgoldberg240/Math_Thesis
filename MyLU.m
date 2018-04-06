% LU Decomposition --- see MATLAB's lu()
% Epperson ch. 7.4, page 430
% Solve AX = I
% Write A = LU -> LUX = I
function x = MyLU(A,p)
%if nargin < 2, [L,U] = LUP(A); end % default to pivoting
if p == 0, x = LUP(A);
else, x = LU(A); end
end
% ------------------- %
function x = LU(a)
% --- Part 1: Factorization --- %
% --- Compute decomposition
n = size(a,1);
for i = 1:n-1
    for j = i+1:n
        a(j,i) = a(j,i)/a(i,i);
        for k = i+1:n
            a(j,k) = a(j,k) - a(j,i)*a(i,k);
        end
        b(j) = b(j) - m*b(i)
    end
end

% % --- Part 2: Solution --- %
% % --- Solve Ly = b
x(1) = b(1);
for i = 2:n
    sum = 0.0;
    for j = 1:i-1
        sum = sum + a(j,i)*x(j);
    end
    x(i) = b(i) - sum;
end
% --- Solve Ux = y
x(n) = x(n)/a(n,n);
for i = n-1:-1:1
    sum = 0.0;
    for j = i+1:n
        sum = sum + a(i,j)*x(j);
    end
    x(i) = (x(i)-sum)/a(i,i);
end
% % ------------------- %
% function [L,U] = LUP(A)

end