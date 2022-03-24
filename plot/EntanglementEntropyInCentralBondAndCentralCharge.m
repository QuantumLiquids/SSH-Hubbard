figure;
Lx = 16:8:40;
S = zeros(size(Lx));

D16=[8000,10000, 12000,14001,16000];
trunc_err16 = [3.37e-06,2.82e-06,2.48e-06, 2.21e-06,1.98e-06]; 
S16=[2.9012856,2.9597793, 3.0059837,3.0438886,3.0746601];
plot(trunc_err16, S16,'o');hold on;


selected_fit_data=1:4;
p = fit(trunc_err16(selected_fit_data)',S16(selected_fit_data)','poly2');
S(1)=p.p3;


D24=8000:2000:16000;
trunc_err24=[3.44e-06,2.86e-06,2.50e-06, 2.23e-06,2.01e-6];
S24=[2.9569130, 3.0166578, 3.0638561, 3.1037447,3.1373982];
plot(trunc_err24, S24,'o');hold on;

selected_fit_data=1:5;
p = fit(trunc_err24(selected_fit_data)',S24(selected_fit_data)','poly2');
S(2)=p.p3;



D32=[8000,9000, 10001,12000, 14000, 16000, 17000,18000];%bond dimension set
trunc_err32=[3.70e-06, 3.28e-06, 3.05e-06, 2.61e-06, 2.32e-06, 2.09e-06, 1.99e-06,1.88e-06];
S32=[2.9556568,2.9858359,3.0139436,3.0604568,3.0989104, 3.1351578,3.1496451,3.1629621];
plot(trunc_err32, S32,'o');hold on;

selected_fit_data=[3:6,7];
p = fit(trunc_err32(selected_fit_data)',S32(selected_fit_data)','poly2');
S(3)=p.p3;



D40=[8000,9000,10000,12000,13000, 14000,15000,16000,17000];%bond dimension set
trunc_err40=[3.47e-6,3.12e-6,2.88e-6,2.50e-06, 2.33e-06,2.20e-06,2.099e-6,1.99e-06, 1.92e-06];
S40=[2.9539752,2.9842074,3.0115933,3.0579595,3.0777987,3.0961124,3.1132416,3.1286635,3.1432339];
plot(trunc_err40, S40,'o');hold on;

selected_fit_data=[3,5:7];
p = fit(trunc_err40(selected_fit_data)',S40(selected_fit_data)','poly2');
S(4)=p.p3;




set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$\varepsilon$','Interpreter','latex');
%ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
ylabel('$S(L_x/2)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 


figure;
plot(Lx(1:4),S(1:4),'o');hold on;


p = fit(log(Lx([1:4])'),S([1:4])','poly1');
fprintf('c=%.5f\n',p.p1*6);
x = Lx(1):0.5:Lx(end);
fl=plot(Lx,p.p1*log(Lx) + p.p2,'-.');
set(gca, 'XScale','log', 'YScale','linear');


set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$L_x$','Interpreter','latex');
%ylabel('$|\langle\Delta_s^\dagger(x)\Delta_s(0)\rangle|$','Interpreter','latex');
ylabel('$S(L_x/2)$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 