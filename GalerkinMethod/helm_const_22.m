function Q = helm_const_22(n,k,a,b,c)
% function Q=helm_const_22(n,k,a,b,c)
%
% Subdiagonal elements A(i,j) with |i-j|=1 
% The matrix for Helmholtz equation
% mit stueckweise konstanter Basisfunktion
% 
%  n:      the number of support points for Gaussian quadrature
%  k:      constant of the equation (real or complex)
%  a,b,c:  Punkte von p (alle Stuetystelle am Rand des Gebietes)
%          bi(x)=1 auf [a,b], =0 sonst
%          bj(x)=1 auf [b,c], =0 sonst 
%          a,b,c complex numbers
% 
%  Example: unit circle as area
%  N: the number of support points on the edge ( = length(p))
%  x = (0:N-1)*(2*pi)/N; p = exp(i*x);
%  A(1,2) = helm_const_22(5,3,p(1),p(2),p(3))

[x,w] = gauss(n);
x = (x+1)/2;
w = w/2;
[y,v] = orthog(n);

[s,t] = meshgrid(x);


st1=k*abs((b-c)*(t+1-s)+(2*b-a-c)*(s-1));
z1=-log(k*abs(b-c+(2*b-a-c)*(s-1)./(t+1-s))).*sing(st1);
Q1=(w*z1)*w'/(2*pi);

t2=repmat(y(:),1,n);
z2=sing(k*abs((b-c)*t2+(2*b-a-c)*(t2.*s-t2)));
Q2=v*(y.*(z2*w'))/(2*pi);

st3=k*abs((b-c)*(t-s.*t+1)+(2*b-a-c)*(s-1));
z3=-log(abs(t-s.*t+1)).*sing(st3).*(1-s); %log(t-s.*t+1)
Q3=(w*z3)*w'/(2*pi);

st4=k*abs((b-a)*s-(c-b)*t+a-b);
z4=(1i/4)*besselh(0,st4+1e-12)+log(st4+1e-12).*sing(st4+1e-12)/(2*pi);
Q4=(w*z4)*w';

Q=abs(b-a)*abs(c-b)*(Q1+Q2+Q3+Q4);