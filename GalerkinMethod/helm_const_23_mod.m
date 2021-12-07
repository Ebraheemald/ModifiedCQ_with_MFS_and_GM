function Q=helm_const_23_mod(n,k,a,b,c,d,Mdist_t)
% k = i*s

% function Q=helm_const_23(n,k,a,b,c,d)
%
% A(i,j) mit |i-j|>1 for Helmholtz equation
% mit stueckweise konstanter Basisfunktion
%
%  n:         the number of support points for the Gaussian quadrature
%  k:         constant of the equation (real or complex)
%  a,b,c,d :  Points of p (all support points at the edge of the area)
%             bi(x)=1 on [a,b], =0 sonst
%             bj(x)=1 on [c,d], =0 sonst 
%             a,b,c,d complex numbers
% 
%  Example: unit circle as area
%  N: the number of support points on the edge ( = length(p))
%  x = (0:N-1)*(2*pi)/N; p = exp(i*x);
%  A(1,3)=helm_const_23(5,3,p(1),p(2),p(3),p(4))


D = min([abs(a-[c d]) abs(b-[c d])]);
if abs(besselh(0,D*k)) < 1e-12
    Q = 0;
else
    [x,w]=gauss(n);
    x=(x+1)/2;
    w=w/2;

    [s,t]=meshgrid(x);


    r=abs((b-a)*s-(d-c)*t+a-c);
    z=(1i/4)*besselh(0,1,k*r+1i*1e-12,1).*exp(1i*k*(r-Mdist_t));

    Q=abs(b-a)*abs(d-c)*(w*z)*w';
end
