addpath('../');
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;
figure_directory='../../figure';

Dset=[12000, 14000, 16000, 17000,18000];%bond dimension set
trunc_err = 1e7*[2.65e-06, 2.32e-06, 2.09e-06, 2.00e-06,1.92e-06];
%D14000 old version(come from 16000) 2.52e-06
extrapolation_poly_degree = 2;
selected_fit_data=[1,2,3,4,5];

D = Dset(1);

FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
SingleParticleCorrelationData = jsondecode(fileread(['../../data/single_particle',FileNamePostfix]));
distance = zeros(1,numel(SingleParticleCorrelationData));
for i=1:numel(SingleParticleCorrelationData)
    FermionSite1 = Site2FermionSite(SingleParticleCorrelationData{i}{1}(1),Ly,Np);
    FermionSite2 = Site2FermionSite(SingleParticleCorrelationData{i}{1}(2),Ly,Np);
    distance(i)=(FermionSite2-FermionSite1)/Ly;
end

SingleParticleCorrelation = zeros(numel(Dset), numel(SingleParticleCorrelationData));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    SingleParticleCorrelationData = jsondecode(fileread(['../../data/single_particle',FileNamePostfix]));
    for i=1:numel(SingleParticleCorrelationData)
        SingleParticleCorrelation(j, i) = SingleParticleCorrelationData{i}{2};
    end
end

[distance, I]=sort(distance); SingleParticleCorrelation=SingleParticleCorrelation(:, I);
h = semilogy(distance,abs(SingleParticleCorrelation),'x');hold on;



SingleParticleCorrelation_ex=zeros(size(distance));
fit_x=trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',SingleParticleCorrelation(selected_fit_data,i),'poly2');
    SingleParticleCorrelation_ex(i)=p.p3;
end

h_ex = semilogy(distance, abs(SingleParticleCorrelation_ex),'o');hold on;

fit_x=[2,3,7,11,15];
%  fit_x=[6,7,10,11,14];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(SingleParticleCorrelation_ex(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = fit_x;
loglog(x,exp(p.p2+p.p1*x),'-.');%fitted line
T=text(10,2.5e-3,['$\xi=',num2str(-1/p.p1, 3),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


% p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
% fprintf('Kc=%.5f\n',-p.p1);
% x = fit_x(1):0.5:fit_x(end);
% fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(5,2e-2,['$K_{c}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);

T2=text(5,4e-1,['$U=',num2str(8),', \lambda = 0.3$']);
set(T2,'Interpreter','latex');set(T2,'Fontsize',24);



l=legend([h; h_ex],'$D=12000$','$14000$','$16000$','$17000$','$18000$','$+\infty$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','SouthWest');




set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$|\langle c^\dagger(0)c(r)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

figure_name_eps = 'single_particle_4x32U8lambda0_3.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['Single particle correlation figure saved as: ', figure_path]);