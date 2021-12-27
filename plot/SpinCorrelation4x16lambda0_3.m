clear;
Lx=16; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;
spin_mode='szsz'; %'szsz', 'spsm','smsp'
begin=3;
endx=14;

Dset=[8000, 10000, 12000, 14000,16000];
trunc_err=  1e7*[3.37e-06,2.82e-06,2.45e-06, 2.18e-06,1.98e-06]; %middle bond;


D = Dset(1);

% FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
%     'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
%     'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
SpinCorrelationData = jsondecode(fileread(['../data/', spin_mode,FileNamePostfix]));
distance = zeros(1,numel(SpinCorrelationData));
for i=1:numel(SpinCorrelationData)
    FermionSite1 = Site2FermionSite(SpinCorrelationData{i}{1}(1),Ly,Np);
    FermionSite2 = Site2FermionSite(SpinCorrelationData{i}{1}(2),Ly,Np);
    distance(i)=(FermionSite2-FermionSite1)/Ly;
end

SpinCorrelation = zeros(numel(Dset), numel(SpinCorrelationData));
for j = 1:numel(Dset)
    D = Dset(j);
%     FileNamePostfix=['begin',num2str(begin),'end',num2str(endx),...
%     'ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),...
%     'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    SpinCorrelationData = jsondecode(fileread(['../data/', spin_mode,FileNamePostfix]));
    for i=1:numel(SpinCorrelationData)
        SpinCorrelation(j, i) = SpinCorrelationData{i}{2} * 2;
    end
end


[distance, I]=sort(distance); SpinCorrelation=SpinCorrelation(:, I);
h = semilogy(distance,abs(SpinCorrelation),'x');hold on;



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex');
ylabel('$|\langle s^x(x) s^x(0)\rangle |$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

