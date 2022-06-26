clear;
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[8000];%bond dimension set
trunc_err = 1e7*[1];
%D14000 old version(come from 16000) 2.52e-06

% extrapolation_poly_degree = 2;
% selected_fit_data=[4,5,6,7]-1;
%change to [2,5,6:8] when data become enough


D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsPSaOY',FileNamePostfix]));
distance=zeros(1,numel(A)/4);
for i=1:numel(A)/4
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(Np+1)/Ly;
end



scsPS=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsPSaOY',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsPSbOY',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsPScOY',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsPSdOY',FileNamePostfix]));
    for i=1:numel(A)
        scsPS(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end



%%======Pyx======%%

PS_data_size = numel(A);
Pyx_data_size = numel(A)/4;
Pyyp_data_size = numel(A)/4;

scsyx=scsPS(:, 1:Pyx_data_size);
scsyx_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
if numel(distance) ~= Pyx_data_size
    error("numel(distance) ~= Pyx_data_size \n");
end
% for i=1:numel(distance)
%     p = fit(fit_x(selected_fit_data)',scsyx(selected_fit_data,i),'poly2');
%     scsyx_ex(i)=p.p3;
% end
h1=loglog(distance, abs(scsyx),'-^');hold on;
%%===== P_yy =======%%
scsyy = scsPS(:, Pyx_data_size+1:2*Pyx_data_size);
scsyy_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
% for i=1:numel(distance)
%     p = fit(fit_x(selected_fit_data)',scsyyp(selected_fit_data,i),'poly2');
%     scsyy_ex(i)=p.p3;
% end
scsyy_ex = scsyy;
h0=loglog(distance, scsyy,'-x');hold on;

fit_x = 2:14;
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end

p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-5:0.5:fit_x(end)+10;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(5,1.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


%%======Pyy'======%%
scsyyp = scsPS(:, 2*Pyx_data_size+1:3*Pyx_data_size);
scsyyp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
if numel(distance) ~= Pyyp_data_size
    error("numel(distance) ~= Pyx_data_size \n");
end
% for i=1:numel(distance)
%     p = fit(fit_x(selected_fit_data)',scsyyp(selected_fit_data,i),'poly2');
%     scsyyp_ex(i)=p.p3;
% end
h2=loglog(distance, -scsyyp,'-x');hold on;


%%======Pyy''======%%
scsyypp = scsPS(:, 3*Pyx_data_size+1:4*Pyx_data_size);
scsyypp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;
% for i=1:numel(distance)
%     p = fit(fit_x(selected_fit_data)',scsyypp(selected_fit_data,i),'poly2');
%     scsyypp_ex(i)=p.p3;
% end
h3=loglog(distance, scsyypp,'-s');hold on;



l=legend([h0,h2,h3,h1],'$\Phi_{yy}(x)$',  '$-\Phi_{yy}^{\prime}(x)$', '$\Phi_{yy}^{\prime\prime}(x)$', '$|\Phi_{yx}(x)|$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
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



