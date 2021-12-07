% The exterior problem for the nonconvex domain (semicircles).
% The CQ error as n (number of steps) increasing,
% where dt = T/N and t = dt*(0:N).

% =============================================================
% =============================================================
% T      : Final time
% N      : Number of time steps to compute the soloution 
% Ne    : Number of time steps to compute the exact soloution 
% dt     : Time step for the solution
% dt_e   : Time step for the exact solution
% lag    : Time lag
% M      : boundary points
% Me     : boundary points for the exact soloution
% X      : Test points
% a2     : Angle of incidence
% a      : a=1 Backward Euler(BE), a=2 the second order backward difference
% formula (BDF2), a=3 Trapezoid rule (TR) for multistep methods

% We have chosen these parameters:
%  M = 1000;   Me = 1500 T = 10;
% a2 =  pi/4; lag = 4;  Ne = 4096
% N = [8, 16, 32, 64, 128, 256, 512, 1024],
% and also, we have choesen these for fner discretization in space:
%  M = 2000;   Me = 3000;

clear all;
fprintf('================\n')
fprintf('The exterior problem for the semicircles\n')
fprintf('================\n')

T=10;                % Final time

% ====== Test pts ========================
t = (0:7)/8; z = 4;
X = z*[(-1-1i)+2*t 1-1i+1i*2*t 1+1i-2*t -1+1i-2*1i*t]; 

% ==== Boundary pts =====================
M =100;
p=semicircles_f(M).';

% ==== Boundary pts for the exact =======
Me=round(1.5*M); pe=semicircles_f(Me).';
pp=8; Ne=2^pp; dt_e=T/Ne; t_e=dt_e*(0:Ne);

% ==== The data =========================
N1 = 2.^(3:(pp-2));  % For different time steps
w1 = [1,5];          % Angular frequency
a2 =  pi/4;          % Angle of incidence
lag = 4;             % Time lag
sgm = .7;
mu = @(t) exp(-(t/sgm).^2);
d_x=@(xx) cos(a2)*real(xx)+sin(a2)*imag(xx);

for i=1:length(w1)
    w=w1(i);
    g=@(xx,t)  sin(w*(t-d_x(xx))) .* mu(t-lag-d_x(xx));
    
    % ========== For the exact ==========
    alpha_e = CQ_solve_sys_mod(Ne,T,pe,g,2);
    u_e = CQ_eval_sys_mod(Ne,T,X,alpha_e,2,pe);

    for a=2:3       % To choose BDF2 and TR for the multistep methods
        err1=zeros(length(N1),1); err2=err1;
        subplot(1,2,a-1)
        
        for j=1:length(N1)
            N=N1(j)                                 % Number of the time steps 
            dt=T/N; t=dt*(0:N);                     % Time steps
            
            % == Without changing variables ==
            alpha1 = CQ_solve_sys(N,T,p,g,a);
            u1 = CQ_eval_sys(N,T,X,alpha1,a,p);      % The solution computed at X 
            err1(j)=max(max(abs(u1(:,:)-u_e(:,1:Ne/N:end))));
            
            % == With changing variables =====                   
            alpha2 = CQ_solve_sys_mod(N,T,p,g,a);
            u2 = CQ_eval_sys_mod(N,T,X,alpha2,a,p);  % The solution computed at X
            err2(j)=max(max(abs(u2(:,:)-u_e(:,1:Ne/N:end))));
        end
        if w==1 && a==2
            BDF2_err1_w1=err1;
            BDF2_err2_w1=err2;
        elseif w==5 && a==2
            BDF2_err1_w5=err1;
            BDF2_err2_w5=err2;
        elseif w==1 && a==3
            TR_err1_w1=err1;
            TR_err2_w1=err2;
        elseif w==5 && a==3
            TR_err1_w5=err1;
            TR_err2_w5=err2;
        end
    end
end

figure;
map = get(gca, 'ColorOrder');
hold on

subplot(1,2,1)
loglog(N1,BDF2_err1_w1,'-','LineWidth', 2); hold on;           % Standard CQ
loglog(N1,BDF2_err1_w5,'x-','Color',map(1,:),'LineWidth', 2);  % Standard CQ
loglog(N1,BDF2_err2_w1,'--','Color',map(2,:),'LineWidth', 2);  % Modified CQ
loglog(N1,BDF2_err2_w5,'x--','Color',map(2,:),'LineWidth', 2); % Modified CQ
h = legend(['standard BDF2 (w=' num2str(w1(1)) ')'], ['standard BDF2 (w=' num2str(w1(2)) ')'],...
        ['modified BDF2 (w=' num2str(w1(1)) ')'],['modified BDF2 (w=' num2str(w1(2)) ')']);
set(h,'Interpreter','Latex','FontSize',12)
h = xlabel('$N$'); set(h,'Interpreter','Latex','FontSize',12)
h = ylabel('error'); set(h,'Interpreter','Latex','FontSize',12) 

subplot(1,2,2)
loglog(N1,TR_err1_w1,'-','LineWidth', 2); hold on;           % Standard CQ
loglog(N1,TR_err1_w5,'x-','Color',map(1,:),'LineWidth', 2);  % Standard CQ
loglog(N1,TR_err2_w1,'--','Color',map(2,:),'LineWidth', 2);  % Modified CQ
loglog(N1,TR_err2_w5,'x--','Color',map(2,:),'LineWidth', 2); % Modified CQ
h = legend(['standard TR (w=' num2str(w1(1)) ')'], ['standard TR (w=' num2str(w1(2)) ')'],...
        ['modified TR (w=' num2str(w1(1)) ')'],['modified TR (w=' num2str(w1(2)) ')']);
set(h,'Interpreter','Latex','FontSize',12)
h = xlabel('$N$'); set(h,'Interpreter','Latex','FontSize',12)
h = ylabel('error'); set(h,'Interpreter','Latex','FontSize',12)
set(gcf, 'Position',  [0, 0, 1000, 700])