clear;
Lx=40; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

begin = 8; 
endx= 32;

Dset=[8000,9000,10000,12000,13000, 14000,15000,17000,16001];%bond dimension set
trunc_err=1e7*[3.47e-6,3.12e-6,2.88e-6,2.50e-06, 2.33e-06,2.20e-06,2.0994e-06, 1.99e-06,2.12e-06];
% grow D17000 trun error = 1.89e-06
extrapolation_poly_degree = 2;
selected_fit_data=[1,4,5,6,7,8];%,1,3,10



D=Dset(1);
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];


SpinCorrelationData = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
distance = zeros(1,numel(SpinCorrelationData));
for i=1:numel(SpinCorrelationData)
    FermionSite1 = Site2FermionSite(SpinCorrelationData{i}{1}(1),Ly,Np);
    FermionSite2 = Site2FermionSite(SpinCorrelationData{i}{1}(2),Ly,Np);
    distance(i)=(FermionSite2-FermionSite1)/Ly;
end

SpinCorrelation = zeros(numel(Dset), numel(SpinCorrelationData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];

    SpinCorrelationDataz = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
    SpinCorrelationDatapm = jsondecode(fileread(['../data/spsm',FileNamePostfix]));
    SpinCorrelationDatamp = jsondecode(fileread(['../data/smsp',FileNamePostfix]));
    for i=1:numel(SpinCorrelationDataz)
        SpinCorrelation(j, i) = SpinCorrelationDataz{i}{2} + 1/2*(SpinCorrelationDatapm{i}{2} + SpinCorrelationDatamp{i}{2});
    end
end


[distance, I]=sort(distance); SpinCorrelation=SpinCorrelation(:, I);


ReducedSpinCorrelation_ex=zeros(size(distance));
fit_x=trunc_err;
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);
for i=1:numel(distance)
    if(distance(i) == Lx/2 -1 )
        p = fit(fit_x(selected_fit_data)',SpinCorrelation(selected_fit_data,i),'poly2');
        ReducedSpinCorrelation_ex(i)= p.p3;
        plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
        
        range=confint(p, 0.95);
        error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
        fprintf("error bar for Spin_correlation_ex at %d = %.6f\n", distance(i), error_bar);
        plot( [0.0, fit_x], [ReducedSpinCorrelation_ex(i), SpinCorrelation(:,i)'], 'o');hold on;
        plot( plot_curve_x, plot_curve_y,'-'); hold on;
    end
end

I=find(distance==Lx/2-1);
fprintf("<S S>(Lx/2-1) = %.6f\n",mean(ReducedSpinCorrelation_ex(I)));



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle s^x(x) s^x(0)\rangle |$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

