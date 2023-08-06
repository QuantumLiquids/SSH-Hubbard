addpath('../');
% figure;
Lx=36; Ly=4;
omega = 5; g = 2.4495; Np = 2; U = 8; Numhole = Lx*Ly/12;

Dset=[8000,10000,12000,14000];%bond dimension set
trunc_err = 1e7*[6.0578e-06,4.8089e-06,3.9459e-06, 3.2e-06];
selected_fit_data=[1,2,3,4];
extrapolation_poly_degree = 2;


Dset = Dset(selected_fit_data);
trunc_err = trunc_err(selected_fit_data);
selected_fit_data = 1:numel(selected_fit_data);


D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
ChargeDensityData = jsondecode(fileread(['../../data/nf',FileNamePostfix]));
distance = zeros(1, size(ChargeDensityData,1));
for i=1:numel(distance)
    FermionSite = Site2FermionSite(ChargeDensityData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end
distance = mean(reshape(distance,4,[]));

ChargeDensity = zeros( numel(Dset), numel(distance) );

for j = 1:numel(Dset)
    D=Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    ChargeDensityData = jsondecode(fileread(['../../data/nf',FileNamePostfix]));
     
    ChargeDensity(j, :) = mean(reshape( transpose(ChargeDensityData(:,2)), 4,[]));
    charge_density_avg = mean(ChargeDensityData(:,2));
    if (charge_density_avg-(1-1/12)) > 1e-8
        error("charge density if not right");
    end
end

% ChargeDensity = (ChargeDensity + ChargeDensity(:,end:-1:1))/2;
h_data=plot(distance + 1, ChargeDensity,'x'); hold on;

ChargeDensity_ex = zeros(1, numel(distance) );


fit_x = trunc_err;
error_bar_set = zeros(1, numel(distance));
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',ChargeDensity(selected_fit_data,i),'poly2');
    % range=confint(p, 0.95);
    % error_bar = (range(2,3) - range(1,3))/2;
    % error_bar_set(i) = error_bar;
    % fprintf("error bar for site %d = %.6f\n", distance(i), error_bar);
    ChargeDensity_ex(i)=p.p3;
end
% fprintf("mean error bar = %.6f\n", mean(error_bar));

% ChargeDensity_ex = (ChargeDensity_ex + ChargeDensity_ex(end:-1:1))/2;
h_ex = plot(distance + 1, ChargeDensity_ex,'-o'); hold on;

% cos_fix_x = Lx/4:3*Lx/4-1;
% % cos_fix_x = Lx/4+1:3*Lx/4-2;
% ChargeDensityLymean = ChargeDensity_ex;
% cos_fix_y = ChargeDensityLymean(cos_fix_x + 1);
% modelfun = @(b,x)(b(2)+ b(3).*cos(pi/2*x+b(1)) );
% mdl = fitnlm(cos_fix_x',cos_fix_y',modelfun,[0.1,1-1/8,0.02]);
% b = mdl.Coefficients.Estimate;
% Acdw = b(3);
% fprintf("A_cdw = %.6f\n",Acdw);
% continous_cos_x = min(cos_fix_x):0.01:max(cos_fix_x);
% plot(continous_cos_x + 1, modelfun(b,continous_cos_x),'-');

% l=legend([h_data',h_ex],'$D=9000$', '$14000$','$16000$','$17000$','extrapolation, $D=\infty$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',20);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'Children'),'markersize',7); % Set line width 1.5 pound
xlabel('$x$','Interpreter','latex');
ylabel('$n(x)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

% set(gca,'XLim',[9,24]);
% set(gca,'Ylim',[0.85,0.94]);
set(gcf,'position',[1000,1000,600,450]);


%{
figure;
% ChargeDensity_ex = (ChargeDensity_ex+ChargeDensity_ex(end:-1:1))/2;

plot(distance(1:end/2), ChargeDensity_ex(1:end/2),'o'); hold on;

distance = distance(1:end/2-2);
ChargeDensity_ex = ChargeDensity_ex(1:end/2-2);

ChargeDensity_ex = ChargeDensity_ex( distance > 2 );
distance = distance( distance > 2 );



set(gca, 'Xlim',[1,Lx/2]);


modelfun = @(b,x)(b(5)+ b(3).*cos(2*b(4).*x+b(1)).*x.^(-b(2)/2) );
mdl = fitnlm(distance',ChargeDensity_ex',modelfun,[ -2.1135
    0.8895
   -0.0357
    0.7866
    0.8757]);

sites = distance;
phi = mdl.Coefficients.Estimate(1);
Kc = mdl.Coefficients.Estimate(2);
deltan = mdl.Coefficients.Estimate(3);
b = mdl.Coefficients.Estimate;
sites = sites(1):0.01:sites(end)+1;
plot(sites, modelfun(b,sites),'-');
l=legend('DMRG data', ['fitting, $K_c= ',num2str(Kc),'$']);
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
