% Position (X) values
x = [0.5, 2.5, 6.5, 10.5, 14.5, 19.7, 20.7];

% Warping data (Y) for each condition
y_control =        [0.004, 0.001, 0.001, 0.002, 0.001, 0.002, 0.000];
y_anneal12h =      [-0.089, -0.054, 0.007, 0.031, 0.020, -0.057, -0.060];
y_strip1h_anneal = [-0.0015, 0.02733, 0.06183, 0.06733, 0.05517, -0.00228, -0.00917];
y_strip2h_anneal = [-0.00125, 0.04025, 0.08625, 0.104, 0.09525, 0.02053, -0.00625];

% Function to offset dataset so that y(1) and y(end) = 0
normalizeWarp = @(y) y - ((y(end) - y(1)) / (x(end) - x(1))) * (x - x(1)) - y(1);

% Normalize each dataset
y_control_norm        = normalizeWarp(y_control);
y_anneal12h_norm      = normalizeWarp(y_anneal12h);
y_strip1h_anneal_norm = normalizeWarp(y_strip1h_anneal);
y_strip2h_anneal_norm = normalizeWarp(y_strip2h_anneal);

% Plot all datasets and fits
figure;
hold on;

datasets = {
    'Control', x, y_control_norm, 'o-', 'b';
    '12h Anneal', x, y_anneal12h_norm, 's-', 'r';
    '1h Strip + 12h Anneal', x, y_strip1h_anneal_norm, 'd-', 'k';
};

for i = 1:size(datasets, 1)
    name = datasets{i,1};
    xi = datasets{i,2};
    yi = datasets{i,3};
    markerStyle = datasets{i,4};
    color = datasets{i,5};
    
    % Plot data
    plot(xi, yi, markerStyle, 'DisplayName', name, 'Color', color);
    
    % Quadratic fit
    p = polyfit(xi, yi, 2);
    x_fit = linspace(min(x), max(x), 200);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, '--', 'Color', color, 'HandleVisibility', 'off'); % Hide duplicate legend
end

% Labels and formatting
xlabel('Position Along DT Strip (mm)');
ylabel('Warp Height (mm, Offset)');
title('Warping of DT Strips After Anneal and Oxide Strip');
legend('Location', 'northwest');
grid on;
