function [a,r]=helm_const_coeff(k,p,g,n)
% function  a=helm_const_coeff(k,p,g,n)
% Coefficient vector for Helmholtz equation with boundary conditions u=g
% mit stueckweise konstanten Basisfunktion
%
%    n:   the number of support points for Gaussian quadrature
%    k:   constant of the equation (real or complex)
%    p:   Support stations on the edge of the area in complex numbers (Nx1)
%    g:   Helmholtz equation u = g on the edge of the area (variable k,x)
%
%  Example: unit circle as area
%  N=32; t=(1:N-1)*(2*pi)/N; p=exp(t*i);
% k = 1+1i;
%  g = @(x) (i/4)*besselh(0,k*abs(x-.5));
%  a = helm_const_coeff(k,p,g)
%
%

if (nargin < 4)
    n = 5;
end

r = l2_rhs(g,p,n);          % compute RHS
A=helm_const_matrix(k,p,n); % compute system matrix
a=A\r;                      % solve system
    