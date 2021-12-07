function alpha = solve_CQ_sys(n,T,x,y,g,a)
% function alpha = solve_CQ_sys(n, T, x, y, g, a)
% returns the coefficient alpha.
% for a given a, the function will choose based on the BE, BDF2, or TR

if a==1
    dlt = @(zt) 1-zt;                % BE
elseif a==2
    dlt = @(zt) 1-zt+.5*(1-zt).^2;   % BDF2
else
    dlt = @(zt) 2*((1-zt)./(1+zt));  % TR
end

K=length(y); M=length(x); dt=T/n; ts=dt*(0:n);

lam = 10^(-12/(2*n+1)); lams = lam.^(0:n);
zeta = exp(2*pi*1i*(0:n)/(n+1));
ss = dlt(lam*zeta.^(-1))/dt;

ghat = zeros(M,n+1); 
for i = 1:M    
    ghat(i,:) = fft(lams.*g(x(i),ts));
end

alpha = zeros(K,n/2+1);
[Y,X] = meshgrid(y,x);
r = abs(X-Y);
for k = 0:n/2    
    V = (1/(2*pi))*besselk(0,ss(k+1)*r);
    alpha(:,k+1) = V\ghat(:,k+1);
end