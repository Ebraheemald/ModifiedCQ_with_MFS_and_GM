function u = CQ_eval_sys_mod(n,T,x,alpha,a,p)
% function u = eval_CQ_sys(n, T, x, y, alpha, a)
% returns the solution.
% for a given a, the function will choose based on the BE, BDF2, or TR.

if a==1
    dlt = @(zt) 1-zt;                % BE
elseif a==2
    dlt = @(zt) 1-zt+.5*(1-zt).^2;   % BDF2
else
    dlt = @(zt) 2*((1-zt)./(1+zt));  % TR
end
M = length(x); dt = T/n;

lam = 10^(-12/(2*n+1)); lams = lam.^(0:n);
zeta = exp(2*pi*1i*(0:n)/(n+1));
ss1=lam*zeta.^(-1);
ss2 = dlt(lam*zeta.^(-1))/dt;

u = zeros(M,n+1);
parfor k = 0:n/2
    u(:,k+1)=helm_const_solx_mod(ss2(k+1)*1i,p,alpha(:,k+1),x,dt,ss1(k+1)); % x is a vector of the test points
end

for k = (n/2+1):n
    u(:,k+1) = conj(u(:,n-k+2));
end

for i = 1:M
    u(i,:) = (1./lams).*ifft(u(i,:));
end