Lx = [16,24,32,40];
Kc = [1.8,1.0225,0.847,0.75];
plot(1./Lx, Kc,'o');hold on;


set(gca, 'Xlim', [0,1/16]);hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$1/Lx$','Interpreter','latex');
ylabel('Kc','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

