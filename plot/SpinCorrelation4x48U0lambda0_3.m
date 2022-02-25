clear;
Lx=48; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 0; Numhole = Lx*Ly/8;

Dset=[10000,11000,12000,13000];%bond dimension set
trunc_err=1e7*[3.64e-06,3.32e-06,3.02e-06,2.78e-06];


extrapolation_poly_degree = 2;
selected_fit_data=1:4;

D = Dset(1);

FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
SpinCorrelationData = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
distance = zeros(1,numel(SpinCorrelationData));
for i=1:numel(SpinCorrelationData)
    FermionSite1 = Site2FermionSite(SpinCorrelationData{i}{1}(1),Ly,Np);
    FermionSite2 = Site2FermionSite(SpinCorrelationData{i}{1}(2),Ly,Np);
    distance(i)=(FermionSite2-FermionSite1)/Ly;
end



SpinCorrelation = zeros(numel(Dset), numel(SpinCorrelationData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    SpinCorrelationDataz = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
    SpinCorrelationDatapm = jsondecode(fileread(['../data/spsm',FileNamePostfix]));
    SpinCorrelationDatamp = jsondecode(fileread(['../data/smsp',FileNamePostfix]));
    for i=1:numel(SpinCorrelationDataz)
        SpinCorrelation(j, i) = SpinCorrelationDataz{i}{2} + 1/2*(SpinCorrelationDatapm{i}{2} + SpinCorrelationDatamp{i}{2});
    end
end

[distance, I]=sort(distance); SpinCorrelation=SpinCorrelation(:, I);
h = semilogy(distance,abs(SpinCorrelation),'x');hold on;



ReducedChargeCorrelation_ex=zeros(size(distance));
fit_x=trunc_err;
for i=1:numel(distance)
    p = fit(fit_x',SpinCorrelation(:,i),'poly2');
    ReducedChargeCorrelation_ex(i)=p.p3;
end

semilogy(distance, abs(ReducedChargeCorrelation_ex),'o');hold on;



fit_x=[4,8,12,16];
%  fit_x=[6,7,10,11,14];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(ReducedChargeCorrelation_ex(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = fit_x;
semilogy(x,exp(p.p2+p.p1*x),'-.');%fitted line
T=text(10,2.5e-3,['$\xi=',num2str(-1/p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);

% 
% p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
% fprintf('Kc=%.5f\n',-p.p1);
% x = fit_x(1):0.5:fit_x(end);
% fl=semilogy(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(8,4.5e-3,['$K_{c}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);


l=legend(h,'$D=10000$', '$11000$','$12000$','$13000$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

set(gca, 'Xlim', [0,16]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle S(x) S(0)\rangle |$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

