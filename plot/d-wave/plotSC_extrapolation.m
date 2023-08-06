clear;
%%% ==== Lx = 16  ====%%%
Lx=16; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;
Dset=[8000,10000, 12000,14000,16000,16001];
trunc_err=  1e7*[3.37e-06,2.82e-06,2.49e-06, 2.21e-06,1.98e-06,2.01e-6]; %middle bond

extrapolation_poly_degree = 2;
selected_fit_data=[1:4,6];

Dset = Dset(selected_fit_data);
trunc_err = trunc_err(selected_fit_data);
D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
index_we_need = find(distance == Lx/2-1);
scsyy=zeros(numel(Dset),1);
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    SC_yy_correlation = 0;
    for i=index_we_need
        SC_yy_correlation = SC_yy_correlation + (A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2});
    end
    scsyy(j) = SC_yy_correlation/Ly;
end


fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);

p = fit(fit_x',scsyy,'poly2');
scsyy_ex=p.p3;
plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
range=confint(p, 0.95);
error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
h1 = plot( [fit_x]/1e7, [ scsyy(:,1)'], 'o');hold on;
l1 = plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;
% set(l1,'MarkerEdgeColor',get(h1,'color'));
% set(l1,'MarkerFaceColor',get(h1,'color'));
set(l1,'Color',get(h1,'color'));
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex));

%%% ==== Lx = 24  ====%%%

Lx=24; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;

begin=4;
endx=20;

Dset=[8000,10000,12000, 14000,16000];
trunc_err=1e7* [3.44e-06,2.86e-06,2.50e-06, 2.23e-06,2.01e-6];%middle bond

extrapolation_poly_degree = 2;

D=Dset(1);
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
index_we_need = find(distance == Lx/2-1);
scsyy=zeros(numel(Dset),1);
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    SC_yy_correlation = 0;
    for i=index_we_need
        SC_yy_correlation = SC_yy_correlation + (A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2});
    end
    scsyy(j) = SC_yy_correlation/Ly;
end


fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);

p = fit(fit_x',scsyy,'poly2');
scsyy_ex=p.p3;
plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
range=confint(p, 0.95);
error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
h2 = plot( [fit_x]/1e7, [ scsyy(:,1)'], 'o');hold on;
l2 = plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;
% set(l2,'MarkerEdgeColor',get(h2,'color'));
% set(l2,'MarkerFaceColor',get(h2,'color'));
set(l2,'Color',get(h2,'color'));
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex));

%%% ==== Lx = 32  ====%%%
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;
Dset=[8000,9000, 10001,12000, 14000, 16000, 17000,18000];%bond dimension set
trunc_err = 1e7*[3.70e-06, 3.28e-06, 3.06e-06, 2.65e-06, 2.32e-06, 2.09e-06, 2.00e-06,1.92e-06];
%D14000 old version(come from 16000) 2.52e-06

extrapolation_poly_degree = 2;
selected_fit_data=[4,5,6:7];


Dset = Dset(selected_fit_data);
trunc_err = trunc_err(selected_fit_data);

D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
index_we_need = find(distance == Lx/2-1);
scsyy=zeros(numel(Dset),1);
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    SC_yy_correlation = 0;
    for i=index_we_need
        SC_yy_correlation = SC_yy_correlation + (A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2});
    end
    scsyy(j) = SC_yy_correlation/Ly;
end


fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);

p = fit(fit_x',scsyy,'poly2');
scsyy_ex=p.p3;
plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
range=confint(p, 0.95);
error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
h3 = plot( [fit_x]/1e7, [ scsyy(:,1)'], 'o');hold on;
l3 = plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;
% set(l3,'MarkerEdgeColor',get(h3,'color'));
% set(l3,'MarkerFaceColor',get(h3,'color'));
set(l3,'Color',get(h3,'color'));
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex));

%%% ==== Lx = 40  ====%%%
Lx=40; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;
begin = 8; 
endx= 32;

Dset=[8000,9000,10000,11000,12000,13000, 14000,15000,16000,17000,18000,16001,14001];%bond dimension set
trunc_err=1e7*[3.47e-6,3.12e-6, 2.88e-6,2.67e-6,2.50e-06, 2.33e-06,2.20e-06,2.0994e-06,1.99e-6, 1.99e-06,1.74e-06,2.12e-06,2.37e-06];
% grow D17000 trun error = 1.89e-06
extrapolation_poly_degree = 2;
selected_fit_data=[5,6,7,8,10];


Dset = Dset(selected_fit_data);
trunc_err = trunc_err(selected_fit_data);

D=Dset(1);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
index_we_need = find(distance == Lx/2-1);
scsyy=zeros(numel(Dset),1);
for j = 1:numel(Dset)
    D = Dset(j);
        FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    SC_yy_correlation = 0;
    for i=index_we_need
        SC_yy_correlation = SC_yy_correlation + (A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2});
    end
    scsyy(j) = SC_yy_correlation/Ly;
end


fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);

p = fit(fit_x',scsyy,'poly2');
scsyy_ex=p.p3;
plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
range=confint(p, 0.95);
error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
h4 = plot( [fit_x]/1e7, [ scsyy(:,1)'], 'o');hold on;
l4 = plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;
% set(l4,'MarkerEdgeColor',get(h4,'color'));
% set(l4,'MarkerFaceColor',get(h4,'color'));
set(l4,'Color',get(h4,'color'));
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex));



l=legend([l1,l2,l3,l4],'$L_x = 16$', '$ 24$',  '$ 32$',  '$40$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'Children'),'markersize',7); % Set line width 1.5 pounds
xlabel('truncation error $\varepsilon$','Interpreter','latex');
ylabel('$\Phi_{yy}(r = \frac{L_x}{2} - 1)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
set(gca,'YLim', [0, inf]);

set(gcf,'position',[1000,1000,500,450]);