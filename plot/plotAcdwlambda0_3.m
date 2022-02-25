Lx = [16,24, 32, 40, 48];
%A_cdw = [0.021293, 0.018159,0.017226, 0.016885, 0.013807];
A_cdw = [0.017277, 0.012838, 0.010395, 0.007962, 0.0];
err = [0.000051,0.000255, 0.000257,  0.000664,0.0];
h1 = errorbar(Lx, A_cdw, err,'o');hold on;
set(gca, 'XScale','log', 'YScale','log')

fit_x = Lx;
p = fit(log(Lx([1:4])'),log(A_cdw([1:4])'),'poly1');
fprintf('Kc=%.5f\n',-p.p1*2);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(Lx,exp(p.p2)*(Lx).^p.p1,'-.');
T=text(30,0.012,['$K_{c}=',num2str(-p.p1*2),'$']);
set(T,'Interpreter','latex');
set(T,'Fontsize',28);

%stop at middle to see if measurement different

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
set(h1,'markersize',10);
xlabel('$L_x$','Interpreter','latex');
ylabel('$A_{cdw}$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
