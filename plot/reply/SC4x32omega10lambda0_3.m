addpath('../');
% figure;
Lx=32; Ly=4;
omega = 10; g = 3.4641016; Np = 3; U = 8; Numhole = Lx*Ly/8;


Dset=[8000,10000,12000,14000,16000];%bond dimension set
trunc_err = 1e7*[1.5249e-06,1.2688e-06, 1.0683e-06,9.1798e-07, 7.8e-07]; %(442,443) bond
extrapolation_poly_degree = 2;
selected_fit_data=[2,3,4];
D=Dset(1);
FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g,8),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
distance=zeros(1,numel(A));
for i=1:numel(A)
    distance(i) = (A{i}{1}(3)-A{i}{1}(1))/(2*Np+1)/Ly;
end

scsyy=zeros(numel(Dset),numel(A));
for j = 1:numel(Dset)
    D = Dset(j);
    FileNamePostfix=['ssh',num2str(Ly),'x',num2str(Lx),'U',num2str(U),'g',num2str(g,8),'omega',num2str(omega),'Np',num2str(Np),'hole',num2str(Numhole),'D',num2str(D),'.json'];
    A = jsondecode(fileread(['../../data/scsyya',FileNamePostfix]));
    B = jsondecode(fileread(['../../data/scsyyb',FileNamePostfix]));
    C = jsondecode(fileread(['../../data/scsyyc',FileNamePostfix]));
    D = jsondecode(fileread(['../../data/scsyyd',FileNamePostfix]));
    for i=1:numel(A)
        scsyy(j,i) = A{i}{2}+B{i}{2}+C{i}{2}+D{i}{2};
    end
end

h=loglog(distance,scsyy,'x');hold on;


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

loglog(distance, scsyy_ex,'o');hold on;
I=find(distance==15);
fprintf("<Delta_yy^dag Delta_yy>(Lx/2-1) = %.6f\n",mean(scsyy_ex(I)));


fit_x=[2,3,6,7,10,11,14];
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(distance==fit_x(i));
    fit_y(i)=mean(scsyy_ex(I));
end


% p = fit((fit_x'),log(abs(fit_y')),'poly1');
% fprintf('correlation length=%.5f\n',-1/p.p1);
% x = fit_x;
% loglog(x,exp(p.p2+p.p1*x),'-.');%fitted line
% T=text(10.2,6e-3,['$\xi=',num2str(-1/p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);
% 
% 
p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(6,8e-3,['$K_{sc}=',num2str(-p.p1, 3),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);



l=legend(h,'$D=8000$', '$10000$', '$12000$', '$14000$', '$16000$');
set(l,'Box','off');set(l,'Interpreter','latex');
set(l,'Fontsize',24);
set(l,'Location','SouthWest');

yticks([1e-5,1e-4,1e-3,1e-2])
yticklabels({'1e-5', '1e-4', '1e-3','1e-2'})

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('$\Phi_{yy}(r)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

%=== dump the figure
figure_directory = '../../figure';
figure_name_eps = 'SComega10.eps';
figure_path = fullfile(figure_directory, figure_name_eps);
saveas(gcf, figure_path, 'epsc');
disp(['the superconductivity correlation figure has been saved at ', figure_path]);


