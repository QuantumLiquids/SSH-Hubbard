clear;
Lx = [16,24, 32, 40];
distance = [7, 11,  15, 19];
SpinCorrelation = [-0.002814, 0.000886,-0.000282, 0.000074];
SpinCorrelation = abs(SpinCorrelation);
err = [0.000026, 0.000001, 0.000011,0.000003];
%15: 6.1e-4
% 9: 1.026e-3
h1 = errorbar(distance, SpinCorrelation, err,'o');hold on;
set(gca, 'XScale','linear', 'YScale','log')

p = fit((distance(1:4)'),log(SpinCorrelation(1:4)'),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = distance(1):0.5:distance(end);
fl=semilogy(distance,exp(p.p2+p.p1*distance),'-.');

T=text(8,2.5e-3,['$\xi=',num2str(-1/p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',28);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(h1,'markersize',10);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle S(0) \cdot S(x)\rangle |$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
