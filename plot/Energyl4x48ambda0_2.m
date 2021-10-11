%weird state
D = [8000,10000,12000];
En = [-433.4199, -433.6052, -433.7367 ];
TrunErr = [ 4.06e-06, 3.10e-06, 2.56e-06 ];%Site  657

plot(TrunErr, En,'-x');hold on;


%half-filled stripe
D = [8000, 10000,12000];
En = [-433.5813, -433.767,-433.870];
TrunErr = [4.39e-06, 3.87e-06, 2.99e-06];

plot(TrunErr, En,'-o');hold on;

set(gca,'Xlim',[0,5e-06]);

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$1/D$','Interpreter','latex');
ylabel('$Energy$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 
