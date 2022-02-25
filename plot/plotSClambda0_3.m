clear;
Lx = [16,24, 32, 40, 48];
distance = [7, 11,  15, 19, 23];
SC = [0.002207, 0.000886,0.000648, 0.000375, 1.1e-4];
err = [0.000053, 0.000005, 0.000035,0.000112, 0];
%15: 6.1e-4
% 9: 1.026e-3
h1 = errorbar(distance, SC, err,'o');hold on;
set(gca, 'XScale','log', 'YScale','log')
p = fit(log(distance(1:4)'),log(SC(1:4)'),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = distance(1):0.5:distance(end);
fl=loglog(distance,exp(p.p2)*distance.^p.p1,'-.');
T=text(8,4.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',28);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(h1,'markersize',10);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('singlet pair correlation','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
