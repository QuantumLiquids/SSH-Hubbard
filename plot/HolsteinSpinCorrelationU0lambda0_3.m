clear;
figure
Lx=32; Ly=4;
omega = 5;
g = 2.4495;
Np=3;
% g = 1.414;
% Np = 2;

U = 1; Numhole = Lx*Ly/8;

Dset=[14000];%bond dimension set

extrapolation_poly_degree = 2;
selected_fit_data=[1];

Dset = Dset(selected_fit_data);
selected_fit_data = 1:numel(selected_fit_data);


D = Dset(1);

FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
SpinCorrelationData = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
distance = zeros(1,numel(SpinCorrelationData));
for i=1:numel(SpinCorrelationData)
    FermionSite1 = (SpinCorrelationData{i}{1}(1))/(Np+1);
    FermionSite2 = (SpinCorrelationData{i}{1}(2))/(Np+1);
    distance(i)=(FermionSite2-FermionSite1)/Ly;
end



SpinCorrelation = zeros(numel(Dset), numel(SpinCorrelationData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['holstein',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    SpinCorrelationDataz = jsondecode(fileread(['../data/szsz',FileNamePostfix]));
    SpinCorrelationDatapm = jsondecode(fileread(['../data/spsm',FileNamePostfix]));
    SpinCorrelationDatamp = jsondecode(fileread(['../data/smsp',FileNamePostfix]));
    for i=1:numel(SpinCorrelationDataz)
%         SpinCorrelation(j, i) =  1/2*(SpinCorrelationDatapm{i}{2} + SpinCorrelationDatamp{i}{2});

         SpinCorrelation(j, i) = SpinCorrelationDataz{i}{2} + 1/2*(SpinCorrelationDatapm{i}{2} + SpinCorrelationDatamp{i}{2});
    end
end

[distance, I]=sort(distance); SpinCorrelation=SpinCorrelation(:, I);
% h = semilogy(distance,abs(SpinCorrelation),'x');hold on;



ReducedChargeCorrelation_ex=SpinCorrelation;

distance2 = mean((reshape(distance,4,[])));
ReducedChargeCorrelation_ex2 = mean((reshape(ReducedChargeCorrelation_ex,4,[])));
h_ex = semilogy(distance2, abs(ReducedChargeCorrelation_ex2),'o');hold on;



fit_x=[4,10,15];
%  fit_x=[6,7,10,11,14];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(ReducedChargeCorrelation_ex(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = min(fit_x)-2:0.5:max(fit_x)+2;
semilogy(x,exp(p.p2+p.p1*x),'-.');%fitted line
T=text(8,2.5e-3,['$\xi=',num2str(-1/p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);

% 
% p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
% fprintf('Kc=%.5f\n',-p.p1);
% x = fit_x(1):0.5:fit_x(end);
% fl=semilogy(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(8,4.5e-3,['$K_{c}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);


% l=legend([h;h_ex],'$D=10000$', '$11000$','$12000$','$13000$','$D=\infty$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');

set(gca, 'Xlim', [0,16]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$|\langle F(r)\rangle |$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


set(gcf,'position',[1000,1000,600,450]);


