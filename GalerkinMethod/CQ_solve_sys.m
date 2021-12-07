function alpha = CQ_solve_sys(n,T,x,g,a)
% function u = eval_CQ_sys(n, T, x, y, alpha, a)
% returns the coefficient alpha.
% For a given a, the function will choose based on the BE, BDF2, or TR.

if a==1
    dlt = @(zt) 1-zt;                % BE
elseif a==2
    dlt = @(zt) 1-zt+.5*(1-zt).^2;   % BDF2
else
    dlt = @(zt) 2*((1-zt)./(1+zt));  % TR
end

M=length(x); dt=T/n; ts=dt*(0:n);

lam = 10^(-12/(2*n+1)); lams = lam.^(0:n);
zeta = exp(2*pi*1i*(0:n)/(n+1));
ss = dlt(lam*zeta.^(-1))/dt;

ghat = zeros(M,n+1);
for k=2:n
    gs = l2_rhs(g,x,ts(k));
    ghat(:,k) = lams(k)*gs;
end

for i=1:M
    ghat(i,:)=fft(ghat(i,:));
end

alpha = zeros(M,n/2+1);
parfor k = 0:n/2
    V = helm_const_matrix(ss(k+1)*1i,x);
    alpha(:,k+1) = V\ghat(:,k+1);
end