%%%%%%%%%%%%%%%%%%%%%%
%  fit sin function and Acdw at finite D first,
%  then extrapolate it to D=inf.
%
%
%
%%%%%%%%
clear;
Lx=24; Ly=4;
omega = 5; g = 2.4495; Np = 3; U = 8; Numhole = Lx*Ly/8;


Dset=[8000, 10000, 12000, 14000,16000];
trunc_err=  1e7*[3.37e-06,2.82e-06,2.45e-06, 2.18e-06,2.01e-06]; %middle bond
% trunc_err=1e7*[7.05e-06, 5.81e-06, 4.95e-06, 4.29e-06];%Site  340


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
        error("charge density is not right");
    end
end

ChargeDensity = (ChargeDensity + ChargeDensity(:,end:-1:1))/2;
% plot(distance + 1, ChargeDensity,'-x'); hold on;

cos_fix_x = Lx/4:3*Lx/4-1;

Acdw_set = zeros(1,numel(Dset));
for j = 1:numel(Dset)
    ChargeDensityLymean = mean(reshape(ChargeDensity(j,:),4,[]));
    cos_fix_y = ChargeDensityLymean(cos_fix_x + 1);
    modelfun = @(b,x)(b(2)+ b(3).*cos(b(4).*x+b(1)) );
    mdl = fitnlm(cos_fix_x',cos_fix_y',modelfun,[0.1,1-1/8,0.02,pi/2]);
    b_fit = mdl.Coefficients.Estimate;
    Acdw_set(j) = b_fit(3);
    fprintf("D = %d, A_cdw = %.6f\n",Dset(j), Acdw_set(j));
    continous_cos_x = min(cos_fix_x):0.01:max(cos_fix_x);
%     plot(continous_cos_x + 1, modelfun(b_fit,continous_cos_x),'-');
end


fit_x=trunc_err;
p = fit(fit_x(1:5)', Acdw_set(1:5)', 'poly2');
Acdw = p.p3;
fprintf("A_cdw = %.6f",Acdw);



plot(trunc_err/1e7, Acdw_set,'o');hold on;
continue_trunc_err = 0:max(trunc_err)/10:max(trunc_err);
% plot(continue_trunc_err/1e7, p.p2 + continue_trunc_err * p.p1,'.-');hold on;
plot(continue_trunc_err/1e7,p.p3 +  p.p2 *continue_trunc_err  + continue_trunc_err.^2 * p.p1,'.-');hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);