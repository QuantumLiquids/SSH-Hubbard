clear;
Lx=32; Ly=4;
omega = 5; 
% g = 2.4495;
 g = 1.414;
Np=2;

U = 0; Numhole = Lx*Ly/8;

Dset=[12000];%bond dimension set

D=Dset(1);
FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];

% ****** yy bond pair ********** %
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end
[order_distance, I] = sort(distance);
order_scsyy = scsyy(I);

h0=loglog(order_distance,(order_scsyy),'-x');hold on;



% ******* On site pair
D=Dset(1);
FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
OnSitePairData = jsondecode(fileread(['../data/onsitesc',FileNamePostfix]));
distance=zeros(1,numel(OnSitePairData));
for i=1:numel(OnSitePairData)
    distance(i) = (OnSitePairData{i}{1}(2)-OnSitePairData{i}{1}(1))/(Np+1)/Ly;
end

scs_onsite=zeros(numel(Dset),numel(OnSitePairData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    OnSitePairData = jsondecode(fileread(['../data/onsitesc',FileNamePostfix]));
    for i=1:numel(OnSitePairData)
        scs_onsite(j,i) = OnSitePairData{i}{2};
    end
end

hs=loglog(distance,scs_onsite,'x');hold on;



fit_x=[1,4,10];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scs_onsite(I));
end

p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-1:0.5:fit_x(end)+6;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(6,3e-2,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);

l=legend([hs,h0],'$\Phi_{s}(x)$','$\Phi_{yy}(x)$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
