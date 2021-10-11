figure;
Lx=48; Ly=4;
omega = 5; g =  0.4472; Np = 1; U = 8; Numhole = Lx*Ly/8;


Dset=[10000,12000,14000,16000,18000];

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

fit_x=1e7*[1.15e-6,9.74e-07, 8.89e-07, 8.26e-07, 7.69e-07];%Site278
for i=1:numel(distance)
    p = fit(fit_x(1:3)',ChargeDensity(1:3,i),'poly2');
    ChargeDensity_ex(i)=p.p3;
end

plot(distance, ChargeDensity_ex,'-o'); hold on;




set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
