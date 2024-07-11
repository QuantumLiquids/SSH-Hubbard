Lx=16; Ly=4;
omega = 5; 
g = 2.4495;
Np=0;

U = 8; Numhole = Lx*Ly/8;


Dset=[10000];

trunc_err = 1e7*[3.70e-06];
selected_fit_data=[1];
extrapolation_poly_degree = 2;
D=Dset(1);
% FileNamePostfix=['mf_ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'hole',num2str(Numhole),'D',num2str(D),'.json'];
FileNamePostfix = '.json';
ChargeDensityData = jsondecode(fileread(['nf.json']));
distance = zeros(1, size(ChargeDensityData,1));
for i=1:numel(distance)
    FermionSite = Site2FermionSite(ChargeDensityData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end

ChargeDensity = zeros( numel(Dset), numel(distance) );

for j = 1:numel(Dset)
    D=Dset(j);
    % FileNamePostfix=['mf_ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    ChargeDensityData = jsondecode(fileread(['nf.json']));
    ChargeDensity(j, :) = transpose(ChargeDensityData(:,2));
    charge_density_avg = mean(ChargeDensityData(:,2));
    if (charge_density_avg-0.875) > 1e-8
        error("charge density if not right");
    end
end


% ChargeDensity = (ChargeDensity + ChargeDensity(:,end:-1:1))/2;


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
    plot(ChargeDensityLymean,'-o'); hold on;
%     plot(continous_cos_x + 1, modelfun(b_fit,continous_cos_x),'-');
end


% fit_x=trunc_err;
% if(extrapolation_poly_degree ==2)
%     p = fit(fit_x(selected_fit_data)', Acdw_set(selected_fit_data)', 'poly2');
%     Acdw = p.p3;
% elseif(extrapolation_poly_degree==3)
%     p = fit(fit_x(selected_fit_data)', Acdw_set(selected_fit_data)', 'poly3');
%     Acdw = p.p4;
% end
% fprintf("A_cdw = %.6f\n",Acdw);
% range=confint(p, 0.95);
% error_bar = (range(2,extrapolation_poly_degree) - range(1,extrapolation_poly_degree))/2;
% fprintf("error bar for A_cdw at %d = %.6f\n", distance(i), error_bar);
% 
% plot(trunc_err/1e7, Acdw_set,'o');hold on;
% continue_trunc_err = 0:max(trunc_err)/10:max(trunc_err);
% % plot(continue_trunc_err/1e7, p.p2 + continue_trunc_err * p.p1,'.-');hold on;
% if(extrapolation_poly_degree==2)
%     plot(continue_trunc_err/1e7,p.p3 +  p.p2 *continue_trunc_err  + continue_trunc_err.^2 * p.p1,'.-');hold on;
% elseif(extrapolation_poly_degree==3)
%     plot(continue_trunc_err/1e7,p.p4 + continue_trunc_err.* (p.p3 +  p.p2 *continue_trunc_err  + continue_trunc_err.^2 * p.p1),'.-');hold on;
% end
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

