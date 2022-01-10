Lx=24; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

begin=4;
endx=20;

Dset=[8000,10000,12000, 14000,16000];
trunc_err=1e7* [3.44e-06,2.79e-06, 2.41e-06, 2.17e-06,2.01e-6];%middle bond


D=Dset(1);
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];

    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

scsyy_ex=zeros(size(distance));
fit_x = trunc_err;%middle bond
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);
for i=1:numel(distance)
    if(distance(i) == Lx/2-1)
        p = fit(fit_x([1:5])',scsyy([1:5],i),'poly2');
        scsyy_ex(i)=p.p3;
        range=confint(p, 0.95);
        error_bar = (range(2,3) - range(1,3))/2;
        fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
        plot( [0.0, fit_x], [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
        plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
%         plot_curve_y =p.p4 + plot_curve_x.*(p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2);
        plot( plot_curve_x, plot_curve_y,'-'); hold on;
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
