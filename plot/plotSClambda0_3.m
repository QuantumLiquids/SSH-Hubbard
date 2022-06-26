clear;
figure;
Lx = [16,24, 32, 40];
distance = [7, 11,  15, 19];
SC = [0.001783, 0.000886, 0.000600, 0.000322];
% SC = [0.002034, 0.000886,0.000626,  0.000504];
err = [0.000022, 0.000005, 0.000050,0.000025];
%15: 6.1e-4
% 9: 1.026e-3
h1 = errorbar(distance(1:4), SC(1:4), err(1:4),'o');hold on;
set(gca, 'XScale','log', 'YScale','log')
p = fit(log(distance([1:2,3])'),log(SC([1:2,3])'),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = distance(1):0.5:distance(end)+2;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(8,4.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',28);

set(gca, 'XTick', [7,11,15,19]);
set(gca,'XTickLabel',{'7','11','15','19'});
% set(gca, 'YLim', [1e-4, 2e-3]);
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(h1,'markersize',10);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$\Phi_{yy}(r)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gca, 'Xlim', [6,20]);
set(gcf,'position',[1000,1000,450,350]);

