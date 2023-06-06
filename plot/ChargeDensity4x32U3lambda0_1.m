clear;
Lx=32; Ly=4;
omega = 5; 
g = 1.414;
Np=2;

U = 3; Numhole = Lx*Ly/8;

Dset=[6000];%bond dimension set
% trunc_err=1e7*[3.64e-06,3.32e-06,3.04e-06,2.80e-06];


D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
distance = zeros(1, size(ChargeDensityData,1));
for i=1:numel(distance)
    FermionSite = Site2FermionSite(ChargeDensityData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end

ChargeDensity = zeros( numel(Dset), numel(distance) );

for j = 1:numel(Dset)
    D=Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
    ChargeDensity(j, :) = transpose(ChargeDensityData(:,2));
    charge_density_avg = mean(ChargeDensityData(:,2));
    if (charge_density_avg-0.875) > 1e-8
        error("charge density if not right");
    end
end

% ChargeDensity = (ChargeDensity + ChargeDensity(:,end:-1:1))/2;
plot(distance + 1, ChargeDensity,'-x'); hold on;
% ChargeDensity_ex = zeros(1, numel(distance) );
% 
% 
% fit_x = trunc_err;
% error_bar_set = zeros(1, numel(distance));
% for i=1:numel(distance)
%     p = fit(fit_x(1:4)',ChargeDensity(1:4,i),'poly1');
%     range=confint(p, 0.95);
%     error_bar = (range(2,2) - range(1,2))/2;
%     error_bar_set(i) = error_bar;
%     fprintf("error bar for site %d = %.6f\n", distance(i), error_bar);
%     ChargeDensity_ex(i)=p.p2;
% end
% fprintf("mean error bar = %.6f\n", mean(error_bar));
% 
% ChargeDensity_ex = (ChargeDensity_ex + ChargeDensity_ex(end:-1:1))/2;
% distancemean=mean(reshape(distance,Ly,[]));
% ChargeDensityLymean = mean(reshape(ChargeDensity_ex,4,[]));
% plot(distancemean + 1, ChargeDensityLymean,'o'); hold on;
% 
% 
% cos_fix_x = Lx/4:3*Lx/4-1;
% 
% cos_fix_y = ChargeDensityLymean(cos_fix_x + 1);
% modelfun = @(b,x)(b(2)+ b(3).*cos(b(4).*x+b(1)) );
% mdl = fitnlm(cos_fix_x',cos_fix_y',modelfun,[0.1,1-1/8,0.02,pi/2]);
% b = mdl.Coefficients.Estimate;
% Acdw = b(3);
% fprintf("A_cdw = %.6f\n",Acdw);
% continous_cos_x = min(cos_fix_x):0.01:max(cos_fix_x);
% plot(continous_cos_x + 1, modelfun(b,continous_cos_x),'-');
% 
% 
% set(gca, 'Xlim',[Lx/4,3*Lx/4-1]+1);
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$n(x)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);




%{
figure;
plot(distance(1:end/2), ChargeDensity_ex(1:end/2),'o'); hold on;

distance = distance(1:end/2);
ChargeDensity_ex = ChargeDensity_ex(1:end/2);

ChargeDensity_ex = ChargeDensity_ex( distance > 2 );
distance = distance( distance > 2 );



set(gca, 'Xlim',[1,Lx/2]);


modelfun = @(b,x)(b(5)+ b(3).*cos(2*b(4).*x+b(1)).*x.^(-b(2)/2) );
mdl = fitnlm(distance',ChargeDensity_ex',modelfun,[1,0.2,1,pi/8,0.91])

sites = distance;
phi = mdl.Coefficients.Estimate(1);
Kc = mdl.Coefficients.Estimate(2);
deltan = mdl.Coefficients.Estimate(3);
b = mdl.Coefficients.Estimate;
sites = sites(1):0.01:sites(end);
plot(sites, modelfun(b,sites),'-');
l=legend('Charge density', ['fitting, $K_c= ',num2str(Kc),'$']);
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
%}

