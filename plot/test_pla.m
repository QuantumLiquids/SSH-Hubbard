clear;
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;


D = 4000;
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
FileNamePostfix=['.json'];
ChargeDensityData = jsondecode(fileread(['../data/nf',FileNamePostfix]));
ChargeDensityData = ChargeDensityData(1:end,:);
disp(mean(ChargeDensityData(:,2)));
% ChargeDensity = (ChargeDensity+ChargeDensity(end:-1:1))/2;
% plot(1:numel(ChargeDensity)/4,ChargeDensity,'-o'); hold on;
distance = zeros(size(ChargeDensityData,1),1);
total_Ly = 2*Np + 4; %Ly = 4
for i=1:numel(distance)
    global_site_num = ChargeDensityData(i, 1);
%     global_y = mod(global_site_num, total_Ly);
%     if(global_y < 2*(Np+1)+1) 
%         y = global_y/ (Np+1);'
%     else
%         y = 3;
%     end
    distance(i) = fix(global_site_num/total_Ly);
end
plot(distance, ChargeDensityData(:,2),'-o'); hold on;

% plot(ChargeDensityData(1:4:end,2),'-o');hold on;
% plot(ChargeDensityData(2:4:end,2),'-o');hold on;
% plot(ChargeDensityData(3:4:end,2),'-o');hold on;
% plot(ChargeDensityData(4:4:end,2),'-o');hold on;



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('Charge Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 