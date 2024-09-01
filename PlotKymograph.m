function PlotKymograph(U,x,T,ui)
% Plots a space-time diagram (Kymograph) of u.

imagesc(x,T,U(:,ui)); 
set(gca,'YDir','normal')
xlabel('$x$','interpreter','latex')
ylabel('$t$','interpreter','latex', 'rotation', 0)
colormap turbo
% shading interp

c = colorbar;
ylabel(c, "u", 'Interpreter', 'latex', 'rotation', 0);
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';

set(gca,'TickLabelInterpreter','latex')
set(gca,'fontsize',24);

end