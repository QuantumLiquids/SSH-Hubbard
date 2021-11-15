Lx = [16,24, 32, 40, 48];
distance = [7, 9,  15, 23];
SC = [0.00142, 0.00079,5e-4,  1.08e-4];
h1 = loglog(distance, SC,'o');hold on;

p = fit(log(distance(1:3)'),log(SC(1:3)'),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
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
