function PlotSolution(dims,U,x,T,ui)

if (dims == 1)
    hold on
    plot(x,U(1,ui),'linewidth',2);
    plot(x,U(end/10,ui),'linewidth',2);
    plot(x,U(end/2,ui),'linewidth',2)
    plot(x,U(end,ui),'linewidth',2)

    h = legend('$t=0$','$t=10\%$','$t=50\%$','$t=100\%$','interpreter','latex');

    set(gca,'fontsize',30);

elseif (dims == 2)
    M = 6;
    Is = [(1:M-1)*round(length(T)/M), size(U,1)];
    m = round(sqrt(numel(ui)));
    clims = [min(U(:,ui),[],'all'),max(U(:,ui),[],'all')];

    for i = Is
        nexttile();
        imagesc(reshape(U(i,ui),m,m));
        axis equal
        axis tight
        title(['$t=', num2str(T(i)),'$'],'interpreter','latex');
        c = colorbar;
        ylabel(c, "H", 'interpreter', 'latex', 'rotation', 0);
        c.TickLabelInterpreter = 'latex';
        caxis(clims);
        xlabel('$x$', 'interpreter', 'latex')
        ylabel('$y$', 'interpreter', 'latex', 'rotation', 0)
        set(gca,'fontsize',24);
    end

    clims = [min(U(end,ui),[],'all'),max(U(end,ui),[],'all')];
    figure
    imagesc(reshape(U(i,ui),m,m));
    axis equal
    axis tight
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    ylabel(c, "H", 'interpreter', 'latex', 'rotation', 0);
    caxis(clims);
    xlabel('$x$', 'interpreter', 'latex')
    ylabel('$y$', 'interpreter', 'latex', 'rotation', 0)

    set(gca,'fontsize',40);

end

end