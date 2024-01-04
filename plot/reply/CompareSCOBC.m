%figure 13(a) blue point in PHYSICAL REVIEW B 95, 125125 (2017)
%"Hybrid-space density matrix renormalization group study of the doped two-dimensional Hubbard model"
%https://journals.aps.org/prb/pdf/10.1103/PhysRevB.95.125125
data = [-0.14641288433382016, 0.6205751693650848
0.9516837481698417, 0.05542881833804428
1.97657393850659, 0.004711785659440073
3.00146412884334, 0.0015945400501635352
3.95314787701318, 0.004079378644037382
4.978038067349928, 0.002495667830903504
9.956076134699853, 0.0008511872297788361
12.957540263543194, 0.0004736471470643617
18.96046852122987, 0.00022865027149876453
23.938506588579802, 0.0001215814119412457
27.96486090775989, 0.00007482257136528207
30.966325036603227, 0.00005071970761189998
32.94289897510981, 0.00004391220378398811
37.99414348462665, 0.00002844839232012267
42.02049780380675, 0.00001839296508760256
44.948755490483165, 0.000011294732959084162
46.99853587115666, 0.000010274856830020704
48.97510980966325, 0.000008059853134334469];
 %the data of D=30000 SU(2)
r = round(data(:,1));
sc_corr = data(:,2)/2;
loglog(r,sc_corr);hold on;

% fit_x=[4,10,15,21,34,43];
% fit_x=[15,21,34,43];
fit_x = r(r>4 & r<12)';
fit_y=zeros(size(fit_x));
for i=1:numel(fit_x)
    I = find(r==fit_x(i));
    fit_y(i)=mean(sc_corr(I));
end

p = fit(log(fit_x'),log(abs(fit_y')),'poly1');
fprintf('Ksc=%.5f\n',-p.p1);
x = fit_x(1):0.5:fit_x(end);
fl=loglog(x,exp(p.p2)*x.^p.p1,'-.');
T=text(10,2.5e-3,['$K_{sc}=',num2str(-p.p1),'$']);
set(T,'Interpreter','latex');set(T,'Fontsize',24);


% p = fit((fit_x'),log(abs(fit_y')),'poly1');
% fprintf('correlation length=%.5f\n',-1/p.p1);
% x = fit_x;
% loglog(x,exp(p.p2+p.p1*x),'-.');%fitted line
% T=text(10,2.5e-3,['$\xi=',num2str(-1/p.p1),'$']);
% set(T,'Interpreter','latex');set(T,'Fontsize',24);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$r$','Interpreter','latex');
ylabel('SC correlation','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 