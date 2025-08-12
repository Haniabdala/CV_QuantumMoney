figure;
hold on;

% Define colors
color_counter = [0.8, 0.355, 0.1];  % for P_counter (eta1 + eta2)
color_homodyne = [0.2, 0.4, 0.7];   % different color for homodyne

% Create dummy lines for legend before plotting
h1 = plot3(NaN, NaN, NaN, 'Color', color_counter, 'LineWidth', 1.5);
h2 = plot3(NaN, NaN, NaN, 'Color', color_homodyne, 'LineWidth', 1.5, 'LineStyle', '--');

for t = 1:length(alphaxp_cell)
    alpha1 = alphaxp_cell{t};
    alpha2 = alpha2xp_cell{t};
    eta1 = eta1xx_cell{t};
    eta2 = eta2xp_cell{t};
    eta5 = homodyne_cell{t};

    % Plot P_counter
    if length(alpha1) == length(alpha2) && length(eta1) == length(eta2)
        plot3(alpha2, alpha1, (0.5 * eta1 + 0.5 * eta2), ...
            'LineWidth', 1.5, 'Color', color_counter);
    else
        warning('Length mismatch (P_counter) at trial %d: skipping.', t);
    end

    % Plot homodyne
    if length(alpha1) == length(alpha2) && length(alpha2) == length(eta5)
        plot3(alpha2, alpha1, eta5, ...
            'LineWidth', 1.5, 'Color', color_homodyne, 'LineStyle', '--');
    else
        warning('Length mismatch (homodyne) at trial %d: skipping.', t);
    end
end

% Set labels and legend before turning hold off
xlabel('Alpha 2');
ylabel('Alpha 1');
zlabel('Probability');

title('3D Comparison: P_{counter} and P_{homodyne}');
legend([h1, h2], {'P_{counter}', 'P_{homodyne}'}, 'Location', 'best');

hold off;
grid on;
view(3);
