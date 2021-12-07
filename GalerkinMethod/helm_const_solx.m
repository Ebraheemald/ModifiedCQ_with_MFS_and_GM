function u=helm_const_solx(k,p,alpha,x,n)
% note k = i*s = delta(lam*zt^(-l))/dt
% ztl = lam*zt^(-l)

% function u = helm_const_solx(n,k,p,alpha,x)
% value of the solution at the point x
% for Helmholtz equation with boundary condition u = g
% with partially constant basic function
% 
%    n:      The number of support points for Gaussian quadrature
%    k:      Constant of the equation (real or complex)
%    p:      Support points at the edge of the area with complex numbers (Nx1)
%    alpha:  Coefficient vector (see helm_const_coeff.m) (Nx1)
%    x:      Points outside the area
%   
%  example: unit circle as area
%  N = 32; t = (0:N-1)*(2*pi)/N; p = exp(i*t);
%  k = 1+1i
%  g = @(x) (i/4)*besselh(0,k*abs(x-.5));
%  alpha = helm_const_coeff(k,p,g);
%  u = helm_const_solx(k,p,alpha,2+2i)
%   

if (nargin < 5)
    n = 5;
end

N = length(p);
m = length(x);

[y,v] = gauss(n);
y = (y+1)/2;
v = v/2;

u = zeros(m,1);

for j = 1:N-1
    s_j = k*abs(x-(p(j+1)-p(j))*y-p(j)); %real(s_j)<0 & imag(s_j)>0
    z_j = (1i/4)*besselh(0,s_j);
    z1 = alpha(j)*abs(p(j+1)-p(j))*v*z_j;
    u = u+z1';
end
sN = k*abs(x-(p(1)-p(N))*y-p(N));
zN = (1i/4)*besselh(0,sN);
z2 = alpha(N)*abs(p(1)-p(N))*v*zN;
u = u+z2';
u = conj(u);
