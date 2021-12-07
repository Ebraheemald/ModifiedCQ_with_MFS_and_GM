function A=helm_const_matrix(k,p,n)
% function A=helm_const_matrix(k,p,n)
%
% matrix for Helmholtz equation with partially constant basis function
% 
%   n: the number of support points for Gaussian quadrature
%   k: constant of the equation (real or complex)
%   p: support points at the edge of the area with complex numbers (Nx1)
%       
%  example: unit circle as area
%  N: the number of support points on the edge ( = length(p))
%  x = (0:N-1)*(2*pi)/N; p = exp(i*x);
%  A = helm_const_matrix(3,p)
%   

if (nargin < 3)
    n = 5;
end

m = length(p);
A = zeros(m);

% these are the diagonal entries (\int_\Gamma_i \int_\Gamma_i)
for i = 1:m-1   
    A(i,i) = helm_const_21(n,k,p(i),p(i+1));   
end
A(m,m) = helm_const_21(n,k,p(m),p(1));

% these are panels touching (distance again zero)
for i = 2:m-1    
    A(i,i-1) = helm_const_22(n,k,p(i-1),p(i),p(i+1));
end
A(m,m-1) = helm_const_22(n,k,p(m-1),p(m),p(1));


% now we look at pairs of panels that are not touching

for i = 3:m-1
    for j = 1:i-2        
        A(i,j) = helm_const_23(n,k,p(i),p(i+1),p(j),p(j+1));       
    end
end

for j = 2:m-2
    A(m,j) = helm_const_23(n,k,p(m),p(1),p(j),p(j+1));
end
A(m,1) = helm_const_22(n,k,p(m),p(1),p(2));

A = A+A.'-diag(diag(A));