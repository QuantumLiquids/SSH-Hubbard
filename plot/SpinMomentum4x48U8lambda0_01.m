figure;
Lx=48; Ly=4;
omega = 5; 
g = 0.4472;
Np=1;

U = 8; Numhole = Lx*Ly/8;

Dset=[10000,12000,14000,16000];%bond dimension set

trunc_err=1e7*[1.15e-6,9.74e-07, 8.89e-07, 8.26e-07];%Site278

selected_fit_data=[1:4];
extrapolation_poly_degree = 2;
D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
SpinMomentData = jsondecode(fileread(['../data/sz',FileNamePostfix]));
distance = zeros(1, size(SpinMomentData,1));
for i=1:numel(distance)
    FermionSite = Site2FermionSite(SpinMomentData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end

SpinMoment = zeros( numel(Dset), numel(distance) );

for j = 1:numel(Dset)
    D=Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    SpinMomentData = jsondecode(fileread(['../data/sz',FileNamePostfix]));
    SpinMoment(j, :) = transpose(SpinMomentData(:,2));
    charge_density_avg = mean(SpinMomentData(:,2));
    if (charge_density_avg) > 1e-8
        error("spin moment not right");
    end

end

% ChargeDensity = (ChargeDensity + ChargeDensity(:,end:-1:1))/2;
h = plot(distance + 1, SpinMoment,'x'); hold on;
SpinMoment_ex = zeros(1, numel(distance) );


fit_x = trunc_err;
error_bar_set = zeros(1, numel(distance));
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',SpinMoment(selected_fit_data,i),'poly1');
%     range=confint(p, 0.95);
%     error_bar = (range(2,3) - range(1,3))/2;
%     error_bar_set(i) = error_bar;
%     fprintf("error bar for site %d = %.6f\n", distance(i), error_bar);
    SpinMoment_ex(i)=p.p2;
end
fprintf("mean error bar = %.6f\n", mean(error_bar));

% ChargeDensity_ex = (ChargeDensity_ex + ChargeDensity_ex(end:-1:1))/2;
plot(distance + 1, SpinMoment_ex,'-o'); hold on;


l=legend(h,'$D=10000$', '$12000$','$14000$','$16000$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$n(x)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gcf,'position',[1000,1000,400,350]);
