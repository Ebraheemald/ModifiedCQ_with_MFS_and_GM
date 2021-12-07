function alpha = solve_CQ_sys_mod(n,T,x,y,g,a)
% function alpha = solve_CQ_sys_mod(n, T, x, y, g, a)
% returns the coefficient alpha size M by N+1.
% based on the BE, BDF2, or TR for a given a.

if a==1
    dlt = @(zt) 1-zt;                % BE
elseif a==2
    dlt = @(zt) 1-zt+.5*(1-zt).^2;   % BDF2
else
    dlt = @(zt) 2*((1-zt)./(1+zt));  % TR
end

K = length(y); M = length(x); dt = T/n; t = dt*(0:n);

lam = 10^(-10/(2*n+1)); lams = lam.^(0:n);
zeta = exp(2*pi*1i*(0:n/2)/(n+1));
ss1=lam*zeta.^(-1);
ss2 = dlt(lam*zeta.^(-1))/dt;

ghat = zeros(M,n+1); 
for i = 1:M        
    ghat(i,:) = fft(lams.*g(x(i),t));
end

alpha = zeros(K,n/2+1);

[Y,X] = meshgrid(y,x);
r = abs(X-Y);
parfor k = 0:n/2    
    V = slp_mod(ss1(k+1),ss2(k+1),r,dt,t);
    alpha(:,k+1) = V\ghat(:,k+1);
end