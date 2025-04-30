A = [899, 1025, 11.5, 180]; % mass
B = [508, 526, 130, 259.08]; % wheel dia

% quadratic regression
p = polyfit(A, B, 2);

% generate fit line
A_fit = linspace(min(A), max(A), 100);
B_fit = polyval(p, A_fit);

% calculate r^2
B_pred = polyval(p, A);
SS_res = sum((B - B_pred).^2);
SS_tot = sum((B - mean(B)).^2);
R_squared = 1 - (SS_res / SS_tot);

% plot data
figure;
scatter(A, B, 'bo', 'filled');
hold on;
plot(A_fit, B_fit, 'r-', 'LineWidth', 2);

% labels/titles
xlabel('Rover Masses, $m$ (kg)', 'interpreter','latex');
ylabel('Wheel Diameter, $W$ (mm)', 'interpreter','latex');
legend('Data points', 'Best-fit curve', 'Location', 'Northwest', 'interpreter','latex');
set(gca, "FontName", 'Times')

% print eq
sprintf('Equation: B = %.4fA^2 + %.4fA + %.4f', p(1), p(2), p(3));

% adding eq to plot
eqn_text = sprintf('$W = %.4f m^2 %+0.4f m %+0.4f$\n$R^2 = %.4f$', p(1), p(2), p(3), R_squared);
text(400, 150, eqn_text, 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'white', 'EdgeColor', 'black');
title('\textbf{Wheel Rover Diameter Versus Wheel Rover Mass}', 'Interpreter', 'latex');

