clear;
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
selected_fit_data = 1:numel(selected_fit_data);


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


scsyy_ex=0;
fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);

if(extrapolation_poly_degree == 2)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,1),'poly2');
    scsyy_ex=p.p3;
    plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
elseif(extrapolation_poly_degree == 3)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,1),'poly3');
    scsyy_ex=p.p4;
    plot_curve_y =p.p4 + plot_curve_x.*(p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2);
elseif(extrapolation_poly_degree == 4)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,1),'poly4');
    scsyy_ex=p.p5;
    plot_curve_y =p.p5 + plot_curve_x.*(p.p4 + plot_curve_x.*(p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2));
end
range=confint(p, 0.95);
error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
% plot( [0.0, fit_x]/1e7, [scsyy_ex, scsyy(:,1)'], 'o');hold on;
plot( [fit_x]/1e7, [ scsyy(:,1)'], 'o');hold on;
plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;


fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex));


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('truncation error','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
