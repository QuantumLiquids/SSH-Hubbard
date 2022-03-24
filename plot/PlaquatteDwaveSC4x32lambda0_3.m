Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[9000, 10001,12000, 14000, 16000];%bond dimension set
trunc_err = 1e7*[ 3.28e-06, 3.06e-06, 2.65e-06, 2.32e-06, 2.09e-06];

extrapolation_poly_degree = 2;
selected_fit_data=[1,4,5];
%change to [2,5,6:8] when data become enough

%%===== P_yy =======%%
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
scsyy_ex=zeros(size(distance));
fit_x = trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyy(selected_fit_data,i),'poly2');
    scsyy_ex(i)=p.p3;
end
distance = mean(transpose(reshape(distance,[],4)));
scsyy_ex = mean(transpose(reshape(scsyy_ex,[],4)));
h0=loglog(distance, scsyy_ex,'-o');hold on;

fit_x=[6,7,10,11];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end

p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-5:0.5:fit_x(end)+10;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);


%%======Pyx======%%

scsPS=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsPSa',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsPSb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsPSc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsPSd',FileNamePostfix]));
    for i=1:numel(A)
        scsPS(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

PS_data_size = numel(A);
Pyx_data_size = numel(A)/4;
Pyyp_data_size = numel(A)/4;

scsyx=scsPS(:, 1:Pyx_data_size);
scsyx_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
if numel(distance) ~= Pyx_data_size
    error("numel(distance) ~= Pyx_data_size \n");
end
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyx(selected_fit_data,i),'poly2');
    scsyx_ex(i)=p.p3;
end
h1=plot(distance, abs(scsyx_ex),'-^');hold on;

%%======Pyy'======%%
scsyyp = scsPS(:, Pyx_data_size+1:2*Pyx_data_size);
scsyyp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
if numel(distance) ~= Pyyp_data_size
    error("numel(distance) ~= Pyx_data_size \n");
end
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyyp(selected_fit_data,i),'poly2');
    scsyyp_ex(i)=p.p3;
end
h2=plot(distance, -scsyyp_ex,'-x');hold on;


%%======Pyy''======%%
scsyypp = scsPS(:, 2*Pyx_data_size+1:3*Pyx_data_size);
scsyypp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyypp(selected_fit_data,i),'poly2');
    scsyypp_ex(i)=p.p3;
end
h3=plot(distance, scsyypp_ex,'-s');hold on;

%%======Pyy'''======%%
scsyyppp = scsPS(:, 3*Pyx_data_size+1:4*Pyx_data_size);
scsyyppp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
for i=1:numel(distance)
    p = fit(fit_x(selected_fit_data)',scsyyppp(selected_fit_data,i),'poly2');
    scsyyppp_ex(i)=p.p3;
end
h4=plot(distance,  -scsyypp_ex,'->');hold on;



l=legend([h0,h2,h3,h4,h1],'$\Phi_{yy}(x)$',  '$-\Phi_{yy}^{\prime}(x)$', '$\Phi_{yy}^{\prime\prime}(x)$', '$-\Phi_{yy}^{\prime\prime\prime}(x)$','$\Phi_{yx}(x)$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',18);
set(l,'Location','SouthWest');


set(gca, 'Xlim', [2,15]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('SC correlations','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


set(gcf,'position',[1000,1000,400,350]);



