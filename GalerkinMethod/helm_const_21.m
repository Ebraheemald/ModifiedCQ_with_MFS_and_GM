function Q = helm_const_21(n,k,a,b)
% function Q = helm_const_21(n,k,a,b)
%
% diagonal elements A (i, i) of the matrix for Helmholtz equation
% mit stueckweise konstanten Basisfunktion
%
%  n:    the number of support points for the Gaussian quadrature
%  k:    constant of the equation (real or complex)
%  a,b:  Points of p (support points at the edge of the area)
%        bi(x)=1 on [a,b], =0 otherwise
%        a,b complex numbers
% 
%  Example: unit circle as area
%  N: the number of support points on the edge (= length (p))
%  x = (0:N-1)*(2*pi)/N; p = exp(i*x);
%  A(1,1)=helm_const_21(5,3,p(1),p(2))

[x,w]=gauss(n);
x=(x+1)/2;
w=w/2;
[y,v]=orthog(n);


% besselh(0,1i*11)/4 approx 1e-6 so we cut off anything bigger than that
L = max(imag(k)*abs(b-a)/11,1);
%L = 1;

x = x/L; w = w/L;

z12 = sing(k*abs((b-a)*y/L)).*(1-y/L);
Q1 = (1/L)*v*z12/pi;

st = k*abs(b-a)*x;

z2 = 2*(1-x).*((1i/4)*besselh(0,st)+log(x*L).*sing(st)/(2*pi));

Q2 = w*z2;

Q = abs(b-a)^2*(Q1+Q2);