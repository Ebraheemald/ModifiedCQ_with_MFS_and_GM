function A = helm_const_matrix_mod(k,p,dt,ztl,n)
% function A = helm_const_matrix_mod(k,p,dt,n)

%note k = i*s = delta(lam*zt^(-l))/dt
%ztl = lam*zt^(-l)
if (nargin < 6)
    n = 5;
end

m = length(p);
A = zeros(m);

%these are the diagonal entries (\int_\Gamma_i \int_\Gamma_i)
for i = 1:m-1   
    A(i,i) = helm_const_21(n,k,p(i),p(i+1));   
end
A(m,m) = helm_const_21(n,k,p(m),p(1));

% these are panels touching (distance again zero)
for i = 2:m-1    
    A(i,i-1) = helm_const_22(n,k,p(i-1),p(i),p(i+1));
end
A(m,m-1) = helm_const_22(n,k,p(m-1),p(m),p(1));


%now we look at pairs of panels that are not touching


for i=3:m-1
    for j=1:i-2
        Mdist = min([abs(p(i)-p(j)), abs(p(i)-p(j+1)), ...
            abs(p(i+1)-p(j)), abs(p(i+1)-p(j+1))]); % approximate distance between panels
        Mdist = floor(Mdist/dt);
        A(i,j) = (ztl^Mdist)*helm_const_23_mod(n,k,p(i),p(i+1),p(j),p(j+1),Mdist*dt);
    end
end

for j=2:m-2
    Mdist = min([abs(p(1)-p(j)),abs(p(1)-p(j+1)), ...
           abs(p(m)-p(j)), abs(p(m)-p(j+1))]);      % approximate distance between panels
    Mdist = floor(Mdist/dt);
    A(m,j)=(ztl^Mdist)*helm_const_23_mod(n,k,p(m),p(1),p(j),p(j+1),Mdist*dt);
end
A(m,1) = helm_const_22(n,k,p(m),p(1),p(2));

A=A+A.'-diag(diag(A));