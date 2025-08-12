load('eta1xx_cell.mat');
load('eta2xp_cell.mat');
load('alpha2xp_cell.mat');
load('alphaxp_cell.mat');
load('homodyne_cell.mat');

figure;
hold on;

% Choose a color (e.g., blue)
lineColor = [0.8, 0.355, 0.1];  % RGB for blue (you can change it)

for t = 1:length(alphaxp_cell)
    alpha1 = alphaxp_cell{t};
    alpha2 = alpha2xp_cell{t};
    eta2 = eta2xp_cell{t};
    eta1 = eta1xx_cell{t};

    % Ensure length matches, and plot using the same color
    if length(alpha1) == length(alpha2) && length(alpha2) == length(eta1)
        plot3(alpha2, alpha1, (0.5 * eta1 + 0.5 * eta2), 'LineWidth', 1.5, 'Color', lineColor);  % z varies with alpha2
    else
        warning('Length mismatch at trial %d: skipping plot.', t);
    end
end

hold off;
grid on;

% Add axis labels
xlabel('Alpha 2');
ylabel('Alpha 1');
zlabel('Probability (eta1)');

% Add a title
title('3D Line Plot of eta1 varying with alpha2');

% Add legend for the entire series
legend('P_{counter}', 'Location', 'best');  % Change 'Line Series' to your preferred label

view(3);  % Set 3D view
