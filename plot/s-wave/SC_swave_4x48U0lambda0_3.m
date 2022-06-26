clear;
figure;
Lx=48; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 0; Numhole = Lx*Ly/8;

Dset=[10000,11000,12000,13000,14000,15000];%bond dimension set
trunc_err=1e7*[3.64e-06,3.32e-06,3.04e-06,2.81e-06,2.56e-06,2.38e-06];


extrapolation_poly_degree = 2;
selected_fit_data=1:4;

% ******* On site pair
D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
OnSitePairData = jsondecode(fileread(['../../data/onsitesc',FileNamePostfix]));
distance=zeros(1,numel(OnSitePairData));
for i=1:numel(OnSitePairData)
    distance(i) = (OnSitePairData{i}{1}(2)- OnSitePairData{i}{1}(1))/(2*Np+1)/Ly;
end

scs_onsite=zeros(numel(Dset),numel(OnSitePairData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    OnSitePairData = jsondecode(fileread(['../../data/onsitesc',FileNamePostfix]));
    for i=1:numel(OnSitePairData)
        scs_onsite(j,i) = OnSitePairData{i}{2};
    end
end


% h2=loglog(distance,scs_onsite,'x');hold on;


scs_ex=zeros(size(distance));
fit_x = trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scs_onsite(selected_fit_data,i),'poly2');
    scs_ex(i)=p.p3;
end

h=loglog(distance, scs_ex,'o');hold on;


fit_x=3:10;
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scs_ex(I));
end



p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-1:0.5:fit_x(end)+6;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(6,3e-2,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);

% l=legend('$D=10000$', '$11000$','$12000$','$13000$', '$\infty$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

scsyy_ex=zeros(size(distance));
fit_x = trunc_err;%middle bond
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly2');
    scsyy_ex(i)=p.p3;
end
h2 = loglog(distance, scsyy_ex,'x');hold on;


set(gca, 'Xlim', [2,16]);
% set(gca, 'Ylim', [2*scs_ex(end),scs_ex(2)]);
set(gca, 'Ylim', [1e-4,scs_ex(2)]);


l=legend([h,h2],'$\Phi_s(r)$', '$\Phi_{yy}(r)$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',28);
set(l,'Location','SouthWest');

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
%ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
ylabel('SC correlation','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gca,'Xlim',[2,16]);
set(gcf,'position',[1000,1000,400,350]);

set(gca, 'XTick', [2,5,10,15]);
set(gca,'XTickLabel',{'2','5','10','15'});