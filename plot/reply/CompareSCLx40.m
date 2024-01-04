% Compare the Superconductivity correlation for lambda = 0 or 0.3. 
% Lx = 32 case.
Lx=40; Ly=4;
U = 8;

Numhole = Lx*Ly/8;

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

loglog(sc_data(:,1), sc_data(:,2),'-x' );hold on;
%==SSHH==%
omega = 5; 
g = 2.4495;
Np=3;

begin = 8;
endx= 32;

Dset=[8000,9000,10000,12000,13000, 14000,15000,16000,17000];%bond dimension set
trunc_err=1e7*[3.47e-6,3.12e-6, 2.88e-6,2.50e-06, 2.33e-06,2.20e-06,2.0994e-06,1.99e-6, 2.00e-06];
% grow D17000 trun error = 1.89e-06
extrapolation_poly_degree = 2;
selected_fit_data=[5,6,7,8,10]-1;
Db=Dset(1);
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    Db = Dset(j);
    FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
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
fit_x = trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly2');
    scsyy_ex(i)=p.p3;
end
MeanDistance = mean(transpose(reshape(distance,[],4)));
MeanScsyy_ex = mean(transpose(reshape(scsyy_ex,[],4)));
h1 = loglog(MeanDistance, MeanScsyy_ex,'o');hold on;


fit_x=[2,5:15];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end



p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = 2:19;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);

% l=legend(h,'$D=8000$','$9000$', '$10000$','$12000$','$13000$','$14000$','$15000$','$16000$', '$17000$', '$18000$', '$16001$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


% 
% l=legend([h1,h2],'$\Phi_{yy}$','$\Phi_s$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');


%%===== Set the figure properties ======%%
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('SC correlation','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
%%===== End of setting figure properties ======%%
