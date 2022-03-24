clear;
Lx=48; Ly=4;
omega = 5; 
% g = 0.4472;
g=1;
% g = 1.414;
% g = 2;
% g = 2.4495;
% g = 2.8284;
% if g<0.0001
%     Np=1;
% elseif g<0.5
%     Np=2;
% elseif g<2.5
%     Np=3;
% else
%     Np=4;
% end
Np=2;
U = 1; Numhole = Lx*Ly/8;


D = 8000;
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
FileNamePostfix=['.json'];
ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
ChargeDensityData = ChargeDensityData(1:end,:);
% ChargeDensity = (reshape(ChargeDensityData(:,2),Ly,[]));
disp(mean(ChargeDensityData(:,2)));
% ChargeDensity = (ChargeDensity+ChargeDensity(end:-1:1))/2;
% plot(1:numel(ChargeDensity)/4,ChargeDensity,'-o'); hold on;
distance = zeros(size(ChargeDensityData,1),1);
for i=1:numel(distance)
    FermionSite = Site2FermionSite(ChargeDensityData(i, 1),Ly,Np);
    distance(i) = fix((FermionSite)/Ly);
end
plot(distance, ChargeDensityData(:,2),'-o'); hold on;



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


% figure;
% 
% 
% D = 6000;
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
% 
% A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
% B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
% C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
% D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
% 
% scsxx=zeros(1,numel(A));
% distance=zeros(1,numel(A));
% 
% for i=1:numel(A)
%     scsxx(1,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
%     distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
% end
% loglog(distance,scsxx,'x');hold on;
% 
% 
% % 
% D = 10000;
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'-2.json'];
% 
% A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
% B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
% C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
% D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
% 
% scsxx=zeros(1,numel(A));
% distance=zeros(1,numel(A));
% 
% for i=1:numel(A)
%     scsxx(1,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
%     distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
% end
% loglog(distance,scsxx,'o');hold on;
% 
% 
% 
% D = 10000;
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'Parallel2.json'];
% 
% A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
% B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
% C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
% D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
% 
% scsxx=zeros(1,numel(A));
% distance=zeros(1,numel(A));
% 
% for i=1:numel(A)
%     scsxx(1,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
%     distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
% end
% plot(distance,scsxx,'^');hold on;
% 
% set(gca,'fontsize',24);
% set(gca,'linewidth',1.5);
% set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% xlabel('$r$','Interpreter','latex');
% ylabel('$Pair Density$','Interpreter','latex');
% set(get(gca,'XLabel'),'FontSize',24); 
% set(get(gca,'YLabel'),'FontSize',24); 
% 
% 
% 
