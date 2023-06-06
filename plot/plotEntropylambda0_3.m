%Lx = 16, D=8000,
S16 = [2.9012856, 2.9597793, 3.0084297];
TruncErr16=1e7*[3.37e-06,2.82e-06,2.46e-06];%Site  228

p = fit(TruncErr16',S16','poly1');
Entropy(1)=p.p2;

S24 = [2.9569625, 3.0157049,  3.0626418]; %8000, 10000, 12000
TruncErr24=1e7* [3.44e-06,2.79e-06, 2.41e-06];%Site  330
p = fit(TruncErr24',S24','poly1');
Entropy(2)=p.p2;

S32 = [2.9556568,  3.0134855, 3.0596893, 3.1025985]; %8000, 10000, 12000, 14000
TruncErr32 = 1e7*[3.70e-06,  2.96e-06, 2.55e-06,2.52e-06 ];%Site  442
p = fit(TruncErr32',S32','poly1');
Entropy(3)=p.p2;

S48 =[ 2.9985328, 3.0462997,  3.0859060, 3.1262995, 3.1522742];%10000,12000, 14000, 16000,18000
TruncErr48=[2.67e-06,  2.34e-06, 2.08e-06,  1.95e-06,  1.76e-06];

p = fit(TruncErr48',S48','poly1');
Entropy(4)=p.p2;

Lx = [16,24, 32, 48];
h1 = plot(log(Lx), Entropy,'o'); hold on;


p = fit(log(Lx(1:end)'),Entropy(1:end)','poly1');
plot( log(Lx), p.p1*log(Lx) + p.p2, '-');hold on;
fprintf("c = %.5f", p.p1*6);


set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
set(h1,'markersize',10);
xlabel('$\ln(L_x)$','Interpreter','latex');
ylabel('$S$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

