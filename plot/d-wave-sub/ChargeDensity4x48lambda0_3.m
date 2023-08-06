
Lx=48; Ly=4;
omega = 5; g = 2.4495; Np = 3; U = 8; Numhole = Lx*Ly/8;


Dset=[8000,10000,12000,14000,16000, 18000, 20000];

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
%     disp(mean(ChargeDensityData(:,2)));
%     ChargeDensity = (ChargeDensity+ChargeDensity(end:-1:1))/2;

    plot(distance+1, ChargeDensity(j,:),'x'); hold on;
end

ChargeDensity_ex = zeros(1, numel(distance) );

%fit_x=1e7*[5.90e-6,4.90e-6,4.19e-06,3.70e-06, 3.35e-06, 3.02e-06, 2.62e-06];%Site  657
fit_x=1e7 * [3.16e-06,2.67e-06,  2.34e-06, 2.08e-06,  1.92e-06,  1.79e-06, 1.47e-06];%middle bond
for i=1:numel(distance)
    p = fit(fit_x(3:6)',ChargeDensity(3:6,i),'poly1');
    ChargeDensity_ex(i)=p.p2;
end

ChargeDensity_ex = (ChargeDensity_ex + ChargeDensity_ex(end:-1:1))/2;
plot(distance+1, ChargeDensity_ex,'o'); hold on;



cos_fix_x = 12:35;
ChargeDensityLymean = mean(reshape(ChargeDensity_ex,4,[]));
cos_fix_y = ChargeDensityLymean(cos_fix_x + 1);
modelfun = @(b,x)(b(2)+ b(3).*cos(b(4).*x+b(1)) );
mdl = fitnlm(cos_fix_x',cos_fix_y',modelfun,[0.1,1-1/8,0.02,pi/2]);
b = mdl.Coefficients.Estimate;
Acdw = b(3);
fprintf("A_cdw = %.6f",Acdw);
continous_cos_x = min(cos_fix_x):0.01:max(cos_fix_x);
plot(continous_cos_x+1, modelfun(b,continous_cos_x),'-');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
