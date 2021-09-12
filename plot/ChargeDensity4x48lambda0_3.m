
Lx=48; Ly=4;
omega = 5; g =  2.4495; Np = 3; U = 8; Numhole = Lx*Ly/8;


Dset=[8000];

D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
distance = zeros(size(ChargeDensityData,1),1);
for i=1:numel(distance)
    FermionSite = Site2FermionSite(ChargeDensityData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end


for j = 1:numel(Dset)
    D=Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
%     ChargeDensityData = ChargeDensityData(1:end,:);
%   ChargeDensity = mean(reshape(ChargeDensityData(:,2),Ly,[]));
    ChargeDensity = ChargeDensityData(:,2);
%     disp(mean(ChargeDensityData(:,2)));
    ChargeDensity = (ChargeDensity+ChargeDensity(end:-1:1))/2;

    plot(distance, ChargeDensity,'-o'); hold on;
end



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
