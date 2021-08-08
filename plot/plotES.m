%plot entanglement spectrum
EEData = jsondecode(fileread(['../data/sv_bond793.json']));

data_size = numel(EEData);
qnN = zeros(1,data_size);
qnSz = zeros(1,numel(EEData));
qnBlockSize = zeros(1,numel(EEData));
maxSV = zeros(1,numel(EEData));
for i = 1:data_size
    qnN(i) = EEData{i}{1}(1);
    qnSz(i) = EEData{i}{1}(2);
    maxSV(i) = EEData{i}{2};
    qnBlockSize(i) = EEData{i}{3};
end
semilogy(qnSz, maxSV,'o');hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$Sz$','Interpreter','latex');
ylabel('Single value','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24); 
set(get(gca,'YLabel'),'FontSize',24); 

    