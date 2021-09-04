Lx=48; Ly=4;
omega = 5; 
g = 1.414;
Np=3;

U = 8; Numhole = Lx*Ly/8;
 
 
 
D = 10000;
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];

A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));

scsxx=zeros(1,numel(A));
distance=zeros(1,numel(A));

for i=1:numel(A)
    scsxx(1,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
semilogy(distance,scsxx,'x');hold on;

scsxx_ex = scsxx;

fit_x=[8,13,18,23];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsxx_ex(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = fit_x;
semilogy(x,exp(p.p2+p.p1*x),'-.');%fitted line
%T=text(7,1.0e-4,['$\xi=',num2str(-1/p.p1),'$']);
%set(T,'Interpreter','latex');set(T,'Fontsize',24);


p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = 5:0.1:23;
fl=semilogy(x,exp(p.p2)*x.^p.p1,'-.');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('Pair Density','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 



