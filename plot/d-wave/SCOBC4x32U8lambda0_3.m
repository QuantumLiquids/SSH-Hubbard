
Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[8000,10000,12000,14000, 14500,16000,18000];%bond dimension set
trunc_err = 1e7*[5.83e-06,5.12e-06,4.7494e-06,4.5936e-06,4.10e-6,4.088e-6,3.7738e-06];
% trunc_err = 1e7*[1/6,1/8,1/10,1/12];
extrapolation_poly_degree = 2;
selected_fit_data=[1,2,3,6,7];

Dset=Dset(selected_fit_data);
trunc_err=trunc_err(selected_fit_data);


% % ******* On site pair ** smaller than pairing on bond 
% D=Dset(1);
% FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
% OnSitePairData = jsondecode(fileread(['../../data/onsitesc',FileNamePostfix]));
% distance=zeros(1,numel(OnSitePairData));
% for i=1:numel(OnSitePairData)
%     distance(i) = (OnSitePairData{i}{1}(2)-OnSitePairData{i}{1}(1))/(2*Np+1)/Ly;
% end
% 
% scs_onsite=zeros(numel(Dset),numel(OnSitePairData));
% for j = 1:numel(Dset)
%     D = Dset(j);
%     FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
%     OnSitePairData = jsondecode(fileread(['../../data/onsitesc',FileNamePostfix]));
%     for i=1:numel(OnSitePairData)
%         scs_onsite(j,i) = OnSitePairData{i}{2};
%     end
% end
% 
% hs=loglog(distance,scs_onsite,'x');hold on;



%%======Pyx======%%
global_Ly = (2*Np+1)*Ly-Np;
D=Dset(1);
FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../../data/scsa',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/global_Ly;
end
distance = distance(1:numel(distance)/4);


scsPS=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../../data/scsa',FileNamePostfix]));
    B = jsondecode(fileread(['../../data/scsb',FileNamePostfix]));
    C = jsondecode(fileread(['../../data/scsc',FileNamePostfix]));
    D = jsondecode(fileread(['../../data/scsd',FileNamePostfix]));
    for i=1:numel(A)
        scsPS(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

PS_data_size = numel(A);
Pyx_data_size = numel(A)/4;
Pyyp_data_size = numel(A)/4;

scsyx=scsPS(:, 1:Pyx_data_size);
scsyx_ex = zeros(1, Pyx_data_size);

fit_x = trunc_err;%middle bond
for i=1:numel(distance)
    if extrapolation_poly_degree==2
        p = fit(fit_x',scsyx(:,i),'poly2');
        scsyx_ex(i)=p.p3;
    elseif extrapolation_poly_degree==1
        p = fit(fit_x',scsyx(:,i),'poly1');
        scsyx_ex(i)=p.p2;
    end
end
h1=loglog(distance, -scsyx_ex,'-^');hold on;


%%======Pyy'======%%
scsyyp = (scsPS(:, Pyx_data_size+1:2*Pyx_data_size)+scsPS(:, 3*Pyx_data_size+1:4*Pyx_data_size))/2;
scsyyp_ex = zeros(1, Pyx_data_size);
for i=1:numel(distance)
    if extrapolation_poly_degree==2
        p = fit(fit_x',scsyyp(:,i),'poly2');
        scsyyp_ex(i)=p.p3;
    elseif extrapolation_poly_degree==1
        p = fit(fit_x',scsyyp(:,i),'poly1');
        scsyyp_ex(i)=p.p2;
    end
end
h2=loglog(distance, scsyyp_ex,'-x');hold on;


fit_x=[2:15];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyyp_ex(I));
end



p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1)-5:0.5:fit_x(end)+10;
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
range=confint(p, 0.8);
fprintf('error bar of Kc = %.12f\n', (range(2,1)-range(1,1))/2);

% T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);





%%======Pyy======%%
scsyypp = scsPS(:, 2*Pyx_data_size+1:3*Pyx_data_size);
scsyypp_ex = zeros(1, Pyx_data_size);
fit_x = trunc_err;%middle bond
for i=1:numel(distance)
    if extrapolation_poly_degree==2
        p = fit(fit_x',scsyypp(:,i),'poly2');
        scsyypp_ex(i)=p.p3;
    elseif extrapolation_poly_degree==1
        p = fit(fit_x',scsyypp(:,i),'poly1');
        scsyypp_ex(i)=p.p2;
    end
end
h3=loglog(distance, scsyypp_ex,'-s');hold on;


% %%======Compare with PBC with the same D=====%%
% D=Dset(1);
% FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
% A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
% distance=zeros(1,numel(A));
% for i=1:numel(A)
%     distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
% end
% 
% scsyy=zeros(numel(Dset),numel(A));
% for j = 1:numel(Dset)
%     D = Dset(j);
%     FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
%     A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
%     B = jsondecode(fileread(['../../data/scsyyb',FileNamePostfix]));
%     C = jsondecode(fileread(['../../data/scsyyc',FileNamePostfix]));
%     D = jsondecode(fileread(['../../data/scsyyd',FileNamePostfix]));
%     for i=1:numel(A)
%         scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
%     end
% end
% 
% h4=loglog(distance,scsyy,'s');hold on;


% l=legend([h1,h2,h3,h4], '$-\Phi_{xy}(x)$',  '$\Phi_{yy}^{\prime}(x)$', '$\Phi_{yy}(x)$','$\Phi_{yy}(x)$ in PBC, for comparison');
l=legend([h1,h2,h3], '$-\Phi_{xy}(r)$, OBC',  '$\Phi_{yy}^{\prime}(r)$', '$\Phi_{yy}(r)$');

set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',28);
set(l,'Location','SouthWest');


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2.5); % Set line width 1.5 pounds
set(get(gca,'Children'),'markersize',9); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
% ylabel('SC correlation','Interpreter','latex');
ylabel('$\Phi(r)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
set(gca,'Xlim',[1,16]);
set(gca,'Ylim',[1e-5,2e-2]);
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2]);
set(gcf,'position',[1000,1000,450,350]);
set(gca, 'XTick', [1,2,5,10,15]);
set(gca,'XTickLabel',{'1','2','5','10','15'});