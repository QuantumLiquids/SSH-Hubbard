clear;
marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;
Lx = [16,24, 32, 40];
distance = [7, 11,  15, 19];
SpinCorrelation = [-0.002821, 0.000886,-0.000321, 0.000092];
SpinCorrelation = abs(SpinCorrelation);
err = [0.000042, 0.000001, 0.000005,0.000003];
%15: 6.1e-4
% 9: 1.026e-3
h1 = errorbar(distance, SpinCorrelation, err,'diamond');hold on;
h1.Color = marker_colors{2};
h1.MarkerSize = 6;
set(gca, 'XScale','linear', 'YScale','log')

T=text(8,4.5e-3,['$\propto e^{-r/\xi_s}$']);
set(T,'Interpreter','latex');set(T,'Fontsize',22);

p = fit((distance(1:4)'),log(SpinCorrelation(1:4)'),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = distance(1)-1:0.5:distance(end)+3;
fl=semilogy(x ,exp(p.p2+p.p1*x ),'--');
fl.Color = marker_colors{6};
% T=text(8,2.5e-3,['$\xi=',num2str(-1/p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',28);

set(gca, 'XTick', [7,11,15,19]);
set(gca,'XTickLabel',{'7','11','15','19'});

set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
% ylabel('$|\langle S(0) \cdot S(x)\rangle |$','Interpreter','latex');
ylabel('$|F(r)|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',20); 
set(get(gca,'YLabel'),'FontSize',20); 
set(gca,'Children',[T, h1, fl]);

% set(gca,'Xlim',[2,16]);
set(gca,'Ylim',[1e-5,1e-2]);
set(gca, 'YTick', [1e-5, 1e-4 1e-3,1e-2]);

set(gcf,'position',[1000,500,400,350]);


