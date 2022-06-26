Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[8000];%bond dimension set



% ****** xx bond pair ********** %
D=Dset(1);
FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../data/scsxxa',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsxx=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['OBCssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../data/scsxxa',FileNamePostfix]));
    B = jsondecode(fileread(['../data/scsxxb',FileNamePostfix]));
    C = jsondecode(fileread(['../data/scsxxc',FileNamePostfix]));
    D = jsondecode(fileread(['../data/scsxxd',FileNamePostfix]));
    for i=1:numel(A)
        scsxx(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

h0=semilogy(distance,scsxx,'o');hold on;



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
set(get(gca,'Children'),'markersize',10); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 



