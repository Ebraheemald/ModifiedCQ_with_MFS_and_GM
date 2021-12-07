function Domain_plot_code(X,x,y)
% function Domain_plot_code(X, x, y)
% plots the domain.

if (nargin < 3)
    n = 1;
else
    n=2;
end

L = size(x,1);

map = get(gca, 'ColorOrder');
map([4,5],:) = map([1,2],:); % swap last two rows.
set(gca, 'ColorOrder',map, 'NextPlot','ReplaceChildren')
%=====================
if n==2
    for l = 1:L
        plot(x(l,:),'linewidth',2); hold on;
        plot(y(l,:),'--','linewidth',2);     
        if (l==1)
            plot(X,'k*'); 
        end
    end
else
    for l = 1:L
        plot(x(l,:),'-','linewidth',2); hold on;  
        if (l==1)
            plot(X,'k*'); 
        end
    end
end
%=====================
if n==2
    h = legend('collocation pts','source pts','$X_j$');
else
    h = legend('panels','$X_j$');
end
set(h,'Interpreter','Latex','FontSize',12)
axis equal;
xlim([-4.5 4.5]); ylim([-4.5 4.5])