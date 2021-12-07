% The exterior problem for the unit circle.
% The CQ error as n (number of steps) increasing,
% where dt = T/N and t = dt*(0:N).

%=============================================================
%=============================================================
% T      : Final time
% N      : Number of time steps to compute the soloution 
% Ne    : Number of time steps to compute the exact soloution 
% dt     : Time step for the solution
% dt_e   : Time step for the exact solution
% lag    : Time lag
% M      : Number of collocation points
% Me     : Number of collocation points for the exact soloution
% K      : Number of source points
% Ke     : Number of source points for the exact soloution
% x      : Collocation points
% y      : Source points
% X      : Test points
% R      : Radius of the circle for the source points
% a2     : Angle of incidence
% a      : a=1 Backward Euler(BE), a=2 the second order backward difference
%formula (BDF2), a=3 Trapezoid rule (TR) for multistep methods

% We have chosen these parameters:
% K = 1000; M = 2000; Ke = 1500 Me = 3000
% T = 10; a2 =  -pi/2; lag = 4; Ne = 4096
% N = [8, 16, 32, 64, 128, 256, 512, 1024]

fprintf('================\n')
fprintf('The exterior problem for the unit circle\n')
fprintf('================\n')

T = 10;                % Final time
R = 0.9;               % Radius of the circle for the source points

% ==== Test pts ==========================
t = (0:7)/8; z = 4;
X = z*[(-1-1i)+2*t 1-1i+1i*2*t 1+1i-2*t -1+1i-2*1i*t]; 

% ==== Collocation and source pts for the probelm ====
K = 1000; M = round(2*K);
ys = @(R, K) R*exp(1i*2*pi*(0:K-1)/K); xs = @(M) exp(1i*2*pi*(0:M-1)/M);
y = ys(R,K); x = xs(M).';

% ==== Collocation and source pts for the exact ====
Ke = 1.5*K; Me = round(2*Ke); 
ye = R*exp(1i*2*pi*(0:Ke-1)/Ke); xe = exp(1i*2*pi*(0:Me-1)/Me).';
p=12; Ne=2^p; dt_e=T/Ne; t_e=dt_e*(0:Ne);

% === The data ===========================
N1 = 2.^(3:(p-2));             % For different time steps
w1=[1 5];                      % Angular frequency
a2 =  -pi/2;                   % Angle of incidence
lag = 4;                       % Time lag
sgm = .7;
mu = @(t) exp(-(t/sgm).^2);
d_x = @(xx) cos(a2)*real(xx)+sin(a2)*imag(xx);

for i=1:length(w1)
    w=w1(i);       
    g=@(xx,t)  sin(w*(t-d_x(xx))) .* mu(t-lag-d_x(xx));
    
    %========== For the exact =================
    alpha_e = solve_CQ_sys_mod(Ne,T,xe,ye,g,2);
    u_e = eval_CQ_sys_mod(Ne,T,X,ye,alpha_e,2); 

    for a=2:3     % To choose BDF2 and TR for the multistep methods
        subplot(1,2,a-1)
        err1=zeros(length(N1),1);  err2=err1;
        for j=1:length(N1)
            N=N1(j)            ;                     % Number of the time steps  
            dt=T/N; t=dt*(0:N);                      % Time steps
            
            % == Without changing variables =====
            alpha1 = solve_CQ_sys(N,T,x,y,g,a);
            u1 = eval_CQ_sys(N,T,X,y,alpha1,a);      % The solution computed at X
            err1(j)=max(max(abs(u1(:,:)-u_e(:,1:Ne/N:end))));
            
            % == With changing variables ========= 
            alpha2 = solve_CQ_sys_mod(N,T,x,y,g,a);        
            u2 = eval_CQ_sys_mod(N,T,X,y,alpha2,a);  % The solution computed at X 
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