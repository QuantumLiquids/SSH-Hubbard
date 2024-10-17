clear;
marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;
figure;
Lx = [16,24, 32, 40];
distance = [7, 11,  15, 19];
SC = [0.001783, 0.000886, 0.000600, 0.000322];
% SC = [0.002034, 0.000886,0.000626,  0.000504];
err = [0.000022, 0.000005, 0.000050,0.000025];
%15: 6.1e-4
% 9: 1.026e-3
h1 = errorbar(distance(1:4), SC(1:4), err(1:4),'>');hold on;
h1.MarkerSize = 6;
h1.Color = marker_colors{4};
set(gca, 'XScale','log', 'YScale','log')
p = fit(log(distance([1:2,3])'),log(SC([1:2,3])'),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = distance(1):0.5:distance(end)+2;
fl=loglog(x,exp(p.p2)*x.^(p.p1-0.03),'--');
fl.Color = marker_colors{6};
range=confint(p, 0.95);
fprintf('error bar of Kc = %.12f\n', (range(2,1)-range(1,1))/2);

% T=text(8,4.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
T=text(8,4.5e-3,['$\propto r^{-K_{sc}}$']);
set(T,'Interpreter','latex');set(T,'Fontsize',22);

set(gca, 'XTick', [7,11,15,19]);
set(gca,'XTickLabel',{'7','11','15','19'});
set(gca, 'YLim', [1e-4, 1e-2]);
set(gca,'fontsize',20);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$\Phi_{yy}(r)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',20); 
set(get(gca,'YLabel'),'FontSize',20); 

set(gca,'Children',[T, h1, fl]);

set(gcf,'position',[1000,500,400,350]);
set(gca,'Xlim',[7,19]);
%set(gca,'Ylim',[1e-5, 1e-1]);

