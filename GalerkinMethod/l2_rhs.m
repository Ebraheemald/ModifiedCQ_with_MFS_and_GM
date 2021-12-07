function r = l2_rhs(g,p,t)
% function r = l2_rhs(g,p,n)
% computation of RHS

if (nargin < 4)
    n = 5;
end

N = length(p);
r = zeros(N,1);

[x,w]=gauss(n);
x=(x+1)/2;
w=w/2;

i=1:N-1;
z=g((p(i+1)-p(i)).*x'+p(i),t);
z1=(w*z')';
r=abs(p(i+1)-p(i)).*z1;


zN=g((p(1)-p(N))*x+p(N),t);
r(N)=abs(p(1)-p(N))*w*zN;