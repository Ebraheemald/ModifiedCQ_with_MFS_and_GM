function u = eval_CQ_sys_mod(n,T,x,y,alpha,a)
% function u = eval_CQ_sys_mod(n, T, x, y, alpha, a)
% returns the solution.
% for a given a, the function will choose based on the BE, BDF2, or TR

if a==1
    dlt = @(zt) 1-zt;                % BE
elseif a==2
    dlt = @(zt) 1-zt+.5*(1-zt).^2;   % BDF2
else
    dlt = @(zt) 2*((1-zt)./(1+zt));  % TR
end

M = length(x); dt = T/n; t = dt*(0:n);

lam = 10^(-10/(2*n+1)); lams = lam.^(0:n);
zeta = exp(2*pi*1i*(0:n/2)/(n+1));
ss1=lam*zeta.^(-1);
ss2 = dlt(lam*zeta.^(-1))/dt;

u = zeros(M,n/2+1);

[Y,X] = meshgrid(y,x);
r = abs(X-Y);
parfor k = 0:n/2
    V = slp_mod(ss1(k+1),ss2(k+1),r,dt,t);    
    u(:,k+1) = V*alpha(:,k+1);
end

for k = (n/2+1):n
   u(:,k+1) = conj(u(:,n-k+2));
end

for i = 1:M
    u(i,:) = (1./lams).*ifft(u(i,:));
end
