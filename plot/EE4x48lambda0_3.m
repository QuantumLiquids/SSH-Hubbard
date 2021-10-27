Lx = 48 ;

load('./EEdata/EE4x48lambda_03D8000.mat', 'S8');

bond_num=numel(S8);
vertical_bond=sort([16-9:28:bond_num,28-9:28:bond_num]);
S8=S8(vertical_bond);
S8=S8(1:2:end)+S8(2:2:end)/2;
x=1:1:Lx-1;

% plot(x, S8, 'o');hold on;

x_points = 1/6*log(Lx/pi*sin(x*pi/Lx));
l1=plot(x_points,S8,'-.x');hold on;


load('./EEdata/EE4x48lambda_03D10000.mat', 'S10');


bond_num=numel(S10);
vertical_bond=sort([16-10:28:bond_num,28-10:28:bond_num]);
S10=S10(vertical_bond);
S10=S10(1:2:end)+S10(2:2:end)/2;
x=1:1:Lx-1;

% plot(x, S10, 'o');hold on;

l2=plot(x_points,S10,'-.x');hold on;


load('./EEdata/EE4x48lambda_03D12000.mat', 'S12');


bond_num=numel(S12);
vertical_bond=sort([16-10:28:bond_num,28-10:28:bond_num]);
S12=S12(vertical_bond);
S12=S12(1:2:end)+S12(2:2:end)/2;
x=1:1:Lx-1;

% plot(x, S12, 'o');hold on;

l3=plot(x_points,S12,'-.x');hold on;

load('./EEdata/EE4x48lambda_03D14000.mat', 'S14');


bond_num=numel(S14);
vertical_bond=sort([16-10:28:bond_num,28-10:28:bond_num]);
S14=S14(vertical_bond);
S14=S14(1:2:end)+S14(2:2:end)/2;
x=1:1:Lx-1;

% plot(x, S14, 'o');hold on;


l4=plot(x_points,S14,'-.x');hold on;




selected_points= abs(x_points-0.18203 )<0.001|abs(x_points-0.338887)<0.001 ...
    | abs(x_points-0.39665)<0.001 | abs(x_points-0.430438)<0.0001 ...
    | abs(x_points-0.448634)<0.00001;
selected_S8 = S8( selected_points );
selected_S10 = S10( selected_points );
selected_S12 = S12( selected_points );
selected_S14 = S14( selected_points );

selected_x_points= x_points(selected_points);
selected_Sdata= [selected_S8;selected_S10;selected_S12;selected_S14];
selected_S_ex=zeros(1, numel(selected_S8));
fit_x=1e7*[5.90e-6,4.90e-6,4.19e-06,3.70e-06];%, 3.33e-06];%Site  657
for i=1:numel(selected_S_ex)
    p = fit(fit_x',selected_Sdata(:,i),'poly1');
    selected_S_ex(i)=p.p2;
end

plot(selected_x_points, selected_S_ex,'o');

p=fit(selected_x_points',selected_S_ex','poly1');
plot(selected_x_points, p.p1*selected_x_points+p.p2,'-');
fprintf('Central Charge=%f\n',p.p1);

for i=1:numel(selected_S_ex)
    p = fit(fit_x',selected_Sdata(:,i),'poly2');
    selected_S_ex(i)=p.p3;
end
plot(selected_x_points, selected_S_ex,'o');
p=fit(selected_x_points',selected_S_ex','poly1');
plot(selected_x_points, p.p1*selected_x_points+p.p2,'-');
fprintf('Central Charge=%f\n',p.p1);




% l=legend('$D=8000$', '$10000$','$12000$','$14000$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
% xlabel('$x$','Interpreter','latex');
xlabel('ranghly $1/6\log x$','Interpreter','latex');
ylabel('$S$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
