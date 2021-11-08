Lx = [24, 32, 40, 48];
distance = [11, 15, 23];
SC = [2.2475e-4, 3.5e-4,5.27e-5];
loglog(distance, SC,'o');hold on;

p = fit(log(distance(2:3)'),log(SC(2:3)'),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(distance,exp(p.p2)*distance.^p.p1,'-.');
T=text(8,4.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('singlet pair correlation','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
