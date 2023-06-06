Lx=32; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;
U = 8; Numhole = Lx*Ly/8;

Db = 8000; %we fix bond dimension as 8000, Dec. 12, 2021




%%%%******** A class *******%%%%
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
scsyy=zeros(1,numel(A));
for i=1:numel(A)
    scsyy(i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
end

loglog(distance,scsyy,'x');hold on;




%%%%******** B class *******%%%%
begin=9;
endx=30;
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
scsyy=zeros(1,numel(A));
for i=1:numel(A)
    scsyy(i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
end

loglog(distance,scsyy,'^');hold on;




%%%%******** D class *******%%%%
begin=7;
endx=30;
FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
    'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
    'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(Db),'.json'];
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end
A = jsondecode(fileread(['../data/scsyya',FileNamePostfix]));
B = jsondecode(fileread(['../data/scsyyb',FileNamePostfix]));
C = jsondecode(fileread(['../data/scsyyc',FileNamePostfix]));
D = jsondecode(fileread(['../data/scsyyd',FileNamePostfix]));
scsyy=zeros(1,numel(A));
for i=1:numel(A)
    scsyy(i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
end

loglog(distance,scsyy,'o');hold on;







y = 10.^(-7:0.5:-1.5);
x = 16*ones(1,numel(y));
plot(x, y,'-.');hold on;


l=legend('$A$ class, $x_0=8$', '$B$ class, $x_0=9$', '$D$ class, $x_0=7$','$x=L_x/2$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('distance $x$','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


