Lx = [16,24, 32,  48];
A_cdw = [0.021398, 0.019870,0.019511, 0.017663];
h1 = loglog(Lx, A_cdw,'o');hold on;

p = fit(log(Lx(1:end)'),log(A_cdw(1:end)'),'poly1');
fprintf('Kc=%.5f\n',-p.p1*2);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(Lx,exp(p.p2)*Lx.^p.p1,'-.');
T=text(30,0.02,['$K_{c}=',num2str(-p.p1*2),'$']);
set(T,'Interpreter','latex');
set(T,'Fontsize',28);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
set(h1,'markersize',10);
xlabel('$L_x$','Interpreter','latex');
ylabel('$A_{cdw}$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
