clear;
marker_colors{1} = [019, 103, 131]/256;
marker_colors{2} = [255,158,002] / 256;
marker_colors{3} = [251,056,071] / 256;
marker_colors{4} = [131,064,028] / 256;
marker_colors{5} = [075,116,178] / 256;
marker_colors{6} = [107,112,092] / 256;
Lx = [16,24, 32, 40];
A_cdw = [0.016924, 0.012861, 0.010620,0.009091];
err = [0.000366,0.000255, 0.000981, 0.000784];
h1 = errorbar(Lx(1:4), A_cdw(1:4), err, 'o'); hold on;
h1.MarkerSize = 6;
h1.Color = marker_colors{1};
set(gca, 'XScale', 'log', 'YScale', 'log')

fit_x = Lx;
p = fit(log(Lx([1:4])'), log(A_cdw([1:4])'), 'poly1');
fprintf('Kc=%.5f\n', -p.p1*2);
x = fit_x(1):0.5:fit_x(end);
fl = loglog(Lx, exp(p.p2) * (Lx).^p.p1, '--');
fl.Color = marker_colors{6};
range = confint(p, 0.95);
fprintf('error bar of Kc = %.12f\n', (range(2,1) - range(1,1)) / 2);

T = text(30, 0.012, ['$\propto r^{-K_{c}}$']);
set(T, 'Interpreter', 'latex');
set(T, 'Fontsize', 22);

set(gca, 'XTick', [16, 24, 32, 40]);
set(gca, 'XTickLabel', {'16', '24', '32', '40'});

set(gca, 'fontsize', 20);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2); % Set line width 1.5 pounds
set(h1, 'markersize', 10);
xlabel('$L_x$', 'Interpreter', 'latex');
ylabel('$A_{\mathrm{cdw}}$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 20); 
set(get(gca, 'YLabel'), 'FontSize', 20); 
set(gca, 'Children', [T, h1, fl]);

set(gca, 'Ylim', [5e-3, 5e-2]);
set(gca, 'Xlim', [16, 40]);

% Add common factor annotation
annotation('textbox', [0.15, 0.85, 0.1, 0.1], 'String', '5 \times', ...
    'Interpreter', 'tex', 'FontSize', 16, 'EdgeColor', 'none');

% Adjust YTick and YTickLabel
set(gca, 'YTick', [5e-3, 5e-2]);
set(gca, 'YTickLabel', {'10^{-3}', '10^{-2}'});

set(gcf, 'position', [1000, 500, 400, 350]);