Lx=48; Ly=4;
omega = 5; 
g = 2.4495;
Np=3;

U = 8; Numhole = Lx*Ly/8;

Dset=[8000,10000,12000,14000,16000,18000];%bond dimension set



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
%fit_x=[1/8,1/10,1/12,1/14];%1/D
fit_x=1e7*[5.90e-6,4.90e-6,4.19e-06,3.70e-06, 3.35e-06, 3.02e-06];%Site  657
plot_curve_x = 0.0:(max(fit_x)/1000):max(fit_x);
for i = 1:numel(distance)
    if(distance(i) == 14)
%         p = fit(fit_x',scsyy(:,i),'poly3');
%         scsyy_ex(i)=p.p4;
%         plot( [0.0, fit_x], [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
%         plot_curve_y = p.p4 + p.p3 .* plot_curve_x + p.p2 .* plot_curve_x.^2 + p.p1 .* plot_curve_x.^3;
%         plot( plot_curve_x, plot_curve_y,'-');
        
        p = fit(fit_x',scsyy(:,i),'poly2');
        scsyy_ex(i)=p.p3;
        plot( [0.0, fit_x], [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
        plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
        plot( plot_curve_x, plot_curve_y,'-'); hold on;
        
        p = fit(fit_x(3:5)',scsyy(3:5,i),'poly2');
        scsyy_ex(i)=p.p3;
        plot( [0.0, fit_x], [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
        plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
        plot( plot_curve_x, plot_curve_y,'-');
        
        p = fit(fit_x(3:6)',scsyy(3:6,i),'poly2');
        scsyy_ex(i)=p.p3;
        plot( [0.0, fit_x], [scsyy_ex(i), scsyy(:,i)'], 'o');hold on;
        plot_curve_y =p.p3 + p.p2 .* plot_curve_x + p.p1 .* plot_curve_x.^2;
        plot( plot_curve_x, plot_curve_y,'-');
        break;
    end
end



% 
% fit_x=[6,7,10,11,14,15,18,19];
% fit_y=zeros(size(fit_x));
% for i=1:numel(fit_x)
%     I = find(distance==fit_x(i));
%     fit_y(i)=mean(scsyy_ex(I));
% end
% 
% 
% p = fit((fit_x'),log(abs(fit_y')),'poly1');
% fprintf('correlation length=%.5f\n',-1/p.p1);
% x = fit_x;
% semilogy(x,exp(p.p2+p.p1*x),'-.');%fitted line
% T=text(10.2,6e-3,['$\xi=',num2str(-1/p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);
% 
% 
% p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
% fprintf('Ksc=%.5f\n',-p.p1);
% x = fit_x(1):0.5:fit_x(end);
% fl=semilogy(x,exp(p.p2)*x.^p.p1,'-.');
% T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);
% 
% 
% 
% l=legend(h,'$D=8000$', '$10000$','$12000$','$14000$');
% set(l,'Box','off');set(l,'Interpreter','latex');
% set(l,'Fontsize',24);
% set(l,'Location','SouthWest');



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('Tuncation error','Interpreter','latex');
ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 





