
Lx=48; Ly=4;
omega = 5; g = 2.4495; Np = 3; U = 8; Numhole = Lx*Ly/8;


Dset=[8000,10000,12000,14000];

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
%     ChargeDensityData = ChargeDensityData(1:end,:);
%   ChargeDensity = mean(reshape(ChargeDensityData(:,2),Ly,[]));
    ChargeDensity(j, :) = transpose(ChargeDensityData(:,2));
%     disp(mean(ChargeDensityData(:,2)));
%     ChargeDensity = (ChargeDensity+ChargeDensity(end:-1:1))/2;

    plot(distance, ChargeDensity(j,:),'-x'); hold on;
end

ChargeDensity_ex = zeros(1, numel(distance) );

fit_x=1e7*[5.90e-6,4.90e-6,4.19e-06,3.70e-06];%Site  657
for i=1:numel(distance)
    p = fit(fit_x(1:end)',ChargeDensity(1:end,i),'poly1');
    ChargeDensity_ex(i)=p.p2;
end

plot(distance, ChargeDensity_ex,'o'); hold on;
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


figure;
ChargeDensity_ex = ChargeDensity_ex+ChargeDensity_ex(end:-1:1);
plot(distance(1:end/2), ChargeDensity_ex(1:end/2),'o'); hold on;


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
