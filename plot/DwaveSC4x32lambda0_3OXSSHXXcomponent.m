
clear;
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[12000];%bond dimension set
trunc_err = 1e7*[1];



D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsxxaOX',FileNamePostfix]));
distance=zeros(1,numel(A)/4);
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(Np+1)/Ly;
end



scsxx=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsxxaOX',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsxxbOX',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsxxcOX',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsxxdOX',FileNamePostfix]));
    for i=1:numel(A)
        scsxx(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

loglog(distance, abs(scsxx),'o');hold on;

fit_x = 2:14;
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsxx(I));
end

p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-5:0.5:fit_x(end)+10;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(5,1.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);



% l=legend([h0,h2,h3,h1],'$|\Phi_{yy}(x)|$',  '$|\Phi_{yy}^{\prime}(x)|$', '$|\Phi_{yy}^{\prime\prime}(x)|$', '$|\Phi_{yx}(x)|$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


set(gca, 'Xlim', [2,15]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('SC correlations','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


set(gcf,'position',[1000,1000,400,350]);



