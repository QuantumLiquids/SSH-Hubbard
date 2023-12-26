bz_edge_color = [038,070,083]/256;
allowed_momentum_line_color=[233,196,107]/256;
% fermi_surface_color=[243,162,097]/256;
fermi_surface_color=[230,111,081]/256;
nodal_line_color = [042,157,142]/256;

mu = 0.5;
% Define the edge length of the square
edgeLength = 2*pi;

% Define the heights of the horizontal dash lines
lineHeights = [-pi, -pi/2, 0, pi/2, pi];

% Define the range and resolution for kx and ky
kx = linspace(-pi, pi, 100);
ky = linspace(-pi, pi, 100);

% Calculate the values of the function cos(kx) + cos(ky)
[X, Y] = meshgrid(kx, ky);
Z = cos(X) + cos(Y);

% Plot the square and horizontal dash lines
figure;
hold on;
axis equal;
axis([-pi pi -pi pi]);
rectangle('Position', [-edgeLength/2, -edgeLength/2, edgeLength, edgeLength], 'EdgeColor', bz_edge_color);
for i = 1:length(lineHeights)
    y = lineHeights(i);
    plot([-pi, pi], [y, y], '--', 'Color', allowed_momentum_line_color);
end

% Plot the nodal lines
plot([-pi, pi], [-pi, pi], '-.', 'Color', nodal_line_color);
plot([-pi, pi], [pi, -pi], '-.', 'Color', nodal_line_color);

% Plot the fermi surface
contourMatrix = contourc(kx, ky, Z, [mu, mu]);% 'LineWidth', 1.5, 'LineColor', 'blue');

% Extract the contour line coordinates
xContour = contourMatrix(1, 2:end);
yContour = contourMatrix(2, 2:end);

% Plot the contour line as a closed shape
fill(xContour, yContour, 'w', 'EdgeColor', fermi_surface_color, 'FaceAlpha', 0.5, 'LineWidth', 1.5);

% Set tick values and labels for the X and Y axes
xticks([-pi, -pi/2, 0, pi/2, pi]);
xticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

% Set labels and title
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',3); % Set line width 1.5 pounds
% set(h1,'markersize',10);
xlabel('$k_x$','Interpreter','latex');
ylabel('$k_y$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
% Set grid
grid off;