% Compare the Superconductivity correlation for lambda = 0 or 0.3. 
% Lx = 32 case.
Lx=32; Ly=4;
U = 8;

Numhole = Lx*Ly/8;

%==Hubbard==%
W = 0;         % Amplitude of CDW chemical potential
t2 = 0;        % NNN hopping
phi = 0;       % twist angle in Y-direction
Kyint = 0;     % Y-direction momentum

Dset = [ 28000];
for i = 1:numel(Dset)
    D = Dset(i);
    FileNamePostfix=['Hubbard',num2str(Ly),'x',num2str(Lx),'t2', num2str(t2),'U', num2str(U),'phi', num2str(phi),  'W', num2str(W),'hole',num2str(Numhole),'D',num2str(D),'T',num2str(Kyint),'.json'];
    [distance, yy_bond_sc_correlation, yy_prime_bond_sc_correlation] = read_hybrid_space_SC(FileNamePostfix, Ly);
    yy_bond_sc_correlation = yy_bond_sc_correlation/16;
    loglog(distance, abs(yy_bond_sc_correlation),'-x'); hold on;
end

fit_x = [5,13,14];

fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(yy_bond_sc_correlation(I));
end


p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = fit_x;
loglog(x,exp(p.p2+p.p1*x),'-.');%fitted line
T=text(7.2,1e-3,['$\xi=',num2str(-1/p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(7,2.5e-4,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


%==Hubbard==%

sc_data = [2, 0.0050251545748088225
3, 0.001938011375271073
4, 0.00266248614403985
5, 0.002271547585601245
6, 0.0011116728815539268
7, 0.0007278953843983146
8, 0.0005586328066431621
9, 0.0004175318936560409
10, 0.0005159928433650851
11, 0.0003755888090680076
12, 0.0003204402976245164
13, 0.00025929437974046724
14, 0.000221221629107045
15, 0.00016534490030002818
16, 0.00014106695349712454
17, 0.00012689610031679235
18, 0.00011414877741680032
19, 0.00010543589908346815
20, 0.00008760496806274342
21, 0.00007278953843983146];

% loglog(sc_data(:,1), sc_data(:,2),'-x' );hold on;

%==SSHH==%
omega = 5; 
g = 2.4495;
Np=3;


Dset=[8000,9000, 10001,12000, 14000, 16000, 17000,18000];%bond dimension set
trunc_err = 1e7*[3.70e-06, 3.28e-06, 3.06e-06, 2.65e-06, 2.32e-06, 2.09e-06, 2.01e-06,1.92e-06];
%D14000 old version(come from 16000) 2.52e-06

extrapolation_poly_degree = 2;
selected_fit_data=[3,5,6:7];
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

% h=loglog(distance,scsyy,'x');hold on;


scsyy_ex=zeros(size(distance));

fit_x = trunc_err;%middle bond
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly2');
    scsyy_ex(i)=p.p3;
%     if(distance(i) == 15)
%         range=confint(p, 0.95);
%         error_bar = (range(2,3) - range(1,3))/2;
%         fprintf("error bar for scsyy_ex at %d = %.6f\n", distance(i), error_bar);
%     end
end

scsyy_ex_mean = mean(reshape(scsyy_ex,[],4),2);

loglog(unique(distance), scsyy_ex_mean,'-o');hold on;
I=find(distance==15);
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex(I)));


fit_x=[5,6,7,10,11,14];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end



p = fit((fit_x'),log(abs(fit_y')),'poly1');
fprintf('correlation length=%.5f\n',-1/p.p1);
x = fit_x;
loglog(x,exp(p.p2+p.p1*x),'-.');%fitted line
T=text(10.2,6e-3,['$\xi=',num2str(-1/p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


% 
l=legend('$\lambda = 0$','$\lambda=0.3$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');




%%===== Set the figure properties ======%%
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('SC correlation','Interpreter','latex');
ylabel('$\Phi_{yy}(r)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);


set(gca, 'XTick', [2,4,8,16]);
set(gca,'XTickLabel',{'2','4','8','16'});
%%===== End of setting figure properties ======%%
figure_directory = '../../figure';
figure_name_eps = 'CompareSC.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['comparison of the superconductivity figure has been saved at ', figure_path]);