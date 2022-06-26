clc;
clear;
close all;
Ly=4;
omega = 5;
g = 2.4495;
Np=3;



U_set = [0,2,3,4,6,8];
Dset = [10000, 8000, 8000, 8000,8000,8000];
Lx_set = [48,48,32,48,48,48];

select = [1,2,3,5,6];
U_set = U_set(select);
Dset = Dset(select);
Lx_set = Lx_set(select);

tile=tiledlayout(1,numel(select),'TileSpacing','Compact','Padding','Compact');
for j = 1:numel(U_set)
    nexttile
    U = U_set(j);
    D=Dset(j);
    Lx = Lx_set(j);
    Numhole = Lx*Ly/8;
    % ****** yy bond pair ********** %
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
    distance=zeros(1,numel(A));
    scsyy=zeros(1,numel(A));
    for i=1:numel(A)
        distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
        scsyy(1,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
    distance = mean(transpose(reshape(distance,[],4)));
    scsyy = mean(transpose(reshape(scsyy,[],4)));
    h=loglog(distance,abs(scsyy),'o');hold on;
    
    
    
    % ******* On site pair ****** %
    OnSitePairData = jsondecode(fileread(['../data/onsitesc',FileNamePostfix]));
    distance=zeros(1,numel(OnSitePairData));
    scs_onsite=zeros(1,numel(OnSitePairData));
    for i=1:numel(OnSitePairData)
        distance(i) = (OnSitePairData{i}{1}(2)-OnSitePairData{i}{1}(1))/(2*Np+1)/Ly;
        scs_onsite(1,i) = OnSitePairData{i}{2};
    end
    h2=loglog(distance,abs(scs_onsite),'x');hold on;
    set(gca,'fontsize',24);
    set(gca,'linewidth',1.5);
    set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
    xlabel('$r$','Interpreter','latex');
    if(j==1)
        ylabel(gca,'SC correlation','Interpreter','latex');
        set(get(gca,'YLabel'),'FontSize',24);

        l=legend([h,h2],'$\Phi_{yy}(r)$', '$|\Phi_{s}(r)|$');
        set(l,'Box','off');set(l,'Interpreter','latex');
        set(l,'Fontsize',28);
        set(l,'Location','SouthWest');
    end
    
    T=text(7,8e-3,['$U = ',num2str(U),'$']);
    set(T,'Interpreter','latex');set(T,'Fontsize',28);
    
    set(get(gca,'XLabel'),'FontSize',24);
    set(gca,'Xlim',[0,16]);
   
    set(gca,'Ylim',[1e-7,1e-1]);
    set(gca, 'XTick', [1,2,5,10,15]);
    set(gca,'XTickLabel',{'1','2','5','10','15'});
    set(gca, 'YTick', [1e-7,1e-5,1e-3,1e-1]);
    set(gcf,'position',[1000,1000,450*4,350]);
    %      if(i > 1)
    %      set(gca, 'YTick', []);
    %      end
end




