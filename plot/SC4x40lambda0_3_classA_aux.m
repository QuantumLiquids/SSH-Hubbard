clear;
Lx=40; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

begin = 8; 
endx= 32;

Dset=[8000,9000,10000,12000,13000, 14000,15000,16000,17000,18000,16001,17001];%bond dimension set
trunc_err=1e7*[3.47e-6,3.12e-6,2.88e-6,2.50e-06, 2.33e-06,2.20e-06,2.0994e-06,1.99e-6, 1.89e-06,1.74e-06,2.12e-06,2.00e-06];
%  wait D=16000, truncation error increase from 1.98e-6 to 2.00e-6
extrapolation_poly_degree = 2;
selected_fit_data=[4,5:6,7,12];

Db=Dset(1);
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    Db = Dset(j);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end



scsyy_ex=zeros(size(distance));
fit_x = trunc_err;
% fit_x = [trunc_err,27];
% A = ones(1,88) * 6.5e-6;
% scsyy = [scsyy; A];
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);
for i=1:numel(distance)
    if(distance(i)==Lx/2-1)
        if(extrapolation_poly_degree == 2)
            p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly2');
            scsyy_ex(i)=p.p3;
            plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
        elseif(extrapolation_poly_degree == 3)
            p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly3');
            scsyy_ex(i)=p.p4;
            plot_curve_y =p.p4 + plot_curve_x.*(p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2);
        elseif(extrapolation_poly_degree == 4)
            p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly4');
            scsyy_ex(i)=p.p5;
            plot_curve_y =p.p5 + plot_curve_x.*(p.p4 + plot_curve_x.*(p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2));
        end
        range=confint(p, 0.95);
        error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
        fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
        plot( [0.0, fit_x]/1e7, [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
        plot( plot_curve_x/1e7, plot_curve_y,'-'); hold on;
    end
end

I=find(distance==Lx/2-1);
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex(I)));





set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('truncation error','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 





