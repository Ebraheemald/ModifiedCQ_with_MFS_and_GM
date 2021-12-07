% The MFS for the exterior problem for the unit circle by choosing sol=1.

% The MFS for the exterior problem for the two ellipses by choosing sol=2.

% The GM for the exterior problem for the semicirles by choosing sol=3.

% sol=4, GM for exterior problem for the semicirles with
% finer discretization is used in space by choosing sol=4.

% For more detail to sol=1, 2, 3, or 4 see Detail_sols.m

sols=1:4; % Choose from 1 to 4
for sol=sols
    if sol==1
        load('1_Disk_R09_z4_K1000_M2000_Ke1500_Me3000_p12_X32')
    elseif sol==2
        load('2_TwoEllipses_R09_z4_K2000_M4000_Ke3000_Me6000_p12_X39')
    elseif sol==3
        load('3_GM_semicircles_z4_M1000_Me1500_pi4_pp12_X32')
    elseif sol==4
        load('4_GM_semicircles_z4_M2000_Me3000_pi4_pp12_X32')
    end
    
    if sol==1
        figure
        Domain_plot_code(X,xs(M),ys(R,K))
    elseif sol==2
        figure
        xs1 = [xs*.5-2; xs*.5+2];
        ys1 = [ys*.5-2;ys*.5+2];           
        Domain_plot_code(X,xs1,ys1);  
    else
        figure
        Domain_plot_code(X,semicircles_f(M))
    end
    figure;          % This figure to plot the convergence
    map = get(gca, 'ColorOrder');
    hold on
    
    subplot(1,2,1)   % 1st subplot to plot the convergence
    
    loglog(N1,BDF2_err1_w1,'-','LineWidth', 2); hold on;           % Standard CQ
    loglog(N1,BDF2_err1_w5,'x-','Color',map(1,:),'LineWidth', 2);  % Standard CQ
    
    loglog(N1,BDF2_err2_w1,'--','Color',map(2,:),'LineWidth', 2);  % Modified CQ
    loglog(N1,BDF2_err2_w5,'x--','Color',map(2,:),'LineWidth', 2); % Modified CQ
    
    h = legend(['standard BDF2 (w=' num2str(w1(1)) ')'], ['standard BDF2 (w=' num2str(w1(2)) ')'],...
        ['modified BDF2 (w=' num2str(w1(1)) ')'],['modified BDF2 (w=' num2str(w1(2)) ')']);
    set(h,'Interpreter','Latex','FontSize',12)
    h = xlabel('$N$'); set(h,'Interpreter','Latex','FontSize',12)
    h = ylabel('error'); set(h,'Interpreter','Latex','FontSize',12)
    if sol==1
        ylim([10^(-6) 10^(1)]);
        legend('Location','southwest')
    elseif sol==2
        ylim([10^(-6) 10^(5)]);
    else
        ylim([10^(-6) 10^(0)]);
        legend('Location','southwest')
    end
          
    subplot(1,2,2)   % 2nd subplot to plot the convergence
    loglog(N1,TR_err1_w1,'-','LineWidth', 2); hold on;           % Standard CQ
    loglog(N1,TR_err1_w5,'x-','Color',map(1,:),'LineWidth', 2);  % Standard CQ
    
    loglog(N1,TR_err2_w1,'--','Color',map(2,:),'LineWidth', 2);  % Modified CQ
    loglog(N1,TR_err2_w5,'x--','Color',map(2,:),'LineWidth', 2); % Modified CQ
    
        
    %=====================
    h = legend(['standard TR (w=' num2str(w1(1)) ')'], ['standard TR (w=' num2str(w1(2)) ')'],...
        ['modified TR (w=' num2str(w1(1)) ')'],['modified TR (w=' num2str(w1(2)) ')']);
    set(h,'Interpreter','Latex','FontSize',12)
    h = xlabel('$N$'); set(h,'Interpreter','Latex','FontSize',12)
    h = ylabel('error'); set(h,'Interpreter','Latex','FontSize',12)
    
    set(gcf, 'Position',  [0, 0, 1000, 700])
    if sol==1
        ylim([10^(-6) 10^(1)]);
        legend('Location','southwest')
    elseif sol==2
        ylim([10^(-6) 10^(5)]);
    else
        ylim([10^(-6) 10^(0)]);
        legend('Location','southwest')
    end    
end