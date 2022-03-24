Lx=48; Ly=4;
omega = 5; 
g = 0.4472;
Np=1;

U = 8; Numhole = Lx*Ly/8;

Dset=[10000,12000,14000,16000,18000];%bond dimension set

D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

% h=semilogy(distance,scsyy,'x');hold on;



scsyy_ex=zeros(size(distance));
fit_x=1e7*[1.15e-6,9.74e-07, 8.89e-07, 8.26e-07, 7.69e-07];%Site278
for i=1:numel(distance)
    p = fit(fit_x(1:end)',scsyy(1:end,i),'poly2');
    scsyy_ex(i)=p.p3;
end

semilogy(distance, scsyy_ex,'o');hold on;


fit_x=[10,18];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = [2,23];
semilogy(x,exp(p.p2+p.p1*x),'-.');%fitted line
% T=text(15.2,1.2e-3,['$\xi=',num2str(-1/p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);


p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
% x = [fit_x(1):0.5:23];
% fl=semilogy(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(15,5.0e-4,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);



% l=legend(h,'$D=10000$', '$12000$','$14000$','$16000$','$18000$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$\Phi_{yy}(x)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

set(gcf,'position',[1000,1000,400,350]);
