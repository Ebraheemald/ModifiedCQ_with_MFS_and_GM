function [x, w] = orthog(n)


j = 1:2*n-1;
alpha = .5*ones(1,2*n-1);
beta = .25./(4-(j-1).^(-2));

anu = zeros(1,2*n);
anu(1) = 1;
tmp = 1/2;
for j = 1:2*n-1
  %tmp = factorial(j)*factorial(j)/factorial(2*j);  
  anu(j+1) = ((-1)^j)*tmp/(j*j+j);
  tmp = tmp*(j+1)*(j+1)/((2*j+1)*(2*j+2));
end

sig = zeros(2*n+1);
ltmp = 2*n;
ltmp = ltmp+1;
sig(2,2:ltmp) = anu;

a = zeros(1,n);b = a;
a(1) = alpha(1)+anu(2)/anu(1);
b(1) = 0;
for k = 3:n+1

  ltmp = 2*n-k+3;
  for l=k:ltmp
    sig(k,l) = sig(k-1,l+1)+(alpha(l-1)-a(k-2))*sig(k-1,l)-...
               b(k-2)*sig(k-2,l)+beta(l-1)*sig(k-1,l-1);
  end
  a(k-1) = alpha(k-1)+sig(k,k+1)/sig(k,k)-sig(k-1,k)/sig(k-1,k-1);
  b(k-1) = sig(k,k)/sig(k-1,k-1);
end


J = diag(a)+diag(sqrt(b(2:end)),1)+diag(sqrt(b(2:end)),-1);
[V,Lam] = eig(J);Lam = diag(Lam);
[Lam,ii] = sort(real(Lam)); V = real(V(:,ii));

x = Lam;
w = V(1,:).*V(1,:);

x = x(:);
w = w(:).';
