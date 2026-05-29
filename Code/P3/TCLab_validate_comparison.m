clear;
clc;
close all;

%% Input
orders_to_compare = [2, 6, 9, 14];

script_dir = fileparts(mfilename('fullpath'));
data_folder = fullfile(script_dir, 'matfiles');
model_folder = fullfile(data_folder, 'model');
figs_folder = fullfile(script_dir, '..', 'figs');

%% Load and Validate Selected Orders
load(fullfile(model_folder, "MSE.mat"))

available_orders = 1:length(MSE);
orders_to_compare = orders_to_compare(ismember(orders_to_compare, available_orders));
if isempty(orders_to_compare)
    error('No valid orders selected. Check `orders_to_compare`.');
end

%% Initialize Storage
num_orders = length(orders_to_compare);
colors = lines(num_orders);

%% Load Data for Each Selected Order
for i = 1:num_orders
    n = orders_to_compare(i);
    filename = fullfile(model_folder, sprintf('simulation_results_n%d.mat', n));
    if ~exist(filename, 'file')
        warning('File not found: %s. Skipping order %d.', filename, n);
        continue;
    end
    load(filename, 't', 'Dy2', 'Dy2_sim', 'error', 'accumulated_variance');
    data{i}.t = t;
    data{i}.Dy2 = Dy2;
    data{i}.Dy2_sim = Dy2_sim;
    data{i}.error = error;
    data{i}.accumulated_variance = accumulated_variance;
    data{i}.order = n;
    data{i}.color = colors(i, :);
end

% Remove any failed loads
data = data(~cellfun(@isempty, data));

if isempty(data)
    error('No valid data loaded. Check file paths and orders.');
end

%% Comparison Validation Data
figure('Units', 'normalized');
theme("light");
set(gcf, 'DefaultAxesFontSize', 24);
set(gcf, 'DefaultTextFontSize', 24);
hold on;
grid on;
plot(data{1}.t, data{1}.Dy2, 'LineWidth', 5, 'Color', 'k', 'DisplayName', 'Validation Data');
for i = 1:num_orders
    plot(data{i}.t, data{i}.Dy2_sim, '-', 'LineWidth', 3, 'Color', data{i}.color, 'DisplayName', sprintf('n=%d', data{i}.order));
end
xlabel('Time [s]');
ylabel('\Delta y [°C]');
title('Validation Data Comparison');
legend('Location', 'best', 'FontSize', 40)
xlim([data{1}.t(1), data{1}.t(end)]);
saveas(gcf, fullfile(figs_folder, 'P3 - Comparison Validation Overlaid.png'));

%% Comparison Errors
figure('Units', 'normalized');
theme("light");
set(gcf, 'DefaultAxesFontSize', 24);
set(gcf, 'DefaultTextFontSize', 24);
hold on;
grid on;

max_errors = zeros(1, num_orders);

for i = 1:num_orders
    plot(data{i}.t, data{i}.error, '-', 'LineWidth', 3, 'Color', data{i}.color, 'DisplayName', sprintf('n=%d', data{i}.order));
    max_errors(i) = max(data{i}.error);
end
xlabel('Time [s]');
ylabel('\Delta y Error [°C]');
title('Error Comparison');
legend('Location', 'best', 'FontSize', 40)
xlim([data{1}.t(1), data{1}.t(end)]);
saveas(gcf, fullfile(figs_folder, 'P3 - Comparison Errors Overlaid.png'));

fprintf('\nMax Error (Cº):\n');
for i = 1:num_orders
    fprintf('Order n=%d: %.2f Cº\n', data{i}.order, max_errors(i));
end

%% Comparison Errors abs
figure('Units', 'normalized');
theme("light");
set(gcf, 'DefaultAxesFontSize', 24);
set(gcf, 'DefaultTextFontSize', 24);
hold on;
grid on;

max_errors = zeros(1, num_orders);

for i = 1:num_orders
    plot(data{i}.t, abs(data{i}.error), '-', 'LineWidth', 3, 'Color', data{i}.color, 'DisplayName', sprintf('n=%d', data{i}.order));
    max_errors(i) = max(data{i}.error);
end
xlabel('Time [s]');
ylabel('Absolute Delta y Error [°C]');
title('Absolute Error Comparison');
legend('Location', 'best', 'FontSize', 40)
xlim([data{1}.t(1), data{1}.t(end)]);
saveas(gcf, fullfile(figs_folder, 'P3 - Comparison Absolute Errors Overlaid.png'));

fprintf('\nAbsolute Max Error (Cº):\n');
for i = 1:num_orders
    fprintf('Order n=%d: %.2f Cº\n', data{i}.order, max_errors(i));
end

%% Comparison normalized error
% Compute normalization factor (max absolute value of actual validation data)
figure('Units', 'normalized');
theme("light");
set(gcf, 'DefaultAxesFontSize', 24);
set(gcf, 'DefaultTextFontSize', 24);
hold on;
grid on;

norm_factor = max(abs(data{1}.Dy2)); % Use first dataset's range for consistency
max_errors = zeros(1, num_orders);

for i = 1:num_orders
    % Calculate percentage error: (|error| / max|Dy2|) * 100
    percentage_error = (abs(data{i}.error) / norm_factor) * 100;
    max_errors(i) = max(percentage_error);
    plot(data{i}.t, percentage_error, '-', 'LineWidth', 3, 'Color', data{i}.color, 'DisplayName', sprintf('n=%d', data{i}.order));

end

xlabel('Time [s]');
ylabel('Normalized Error [%]');
title('Normalized Error Comparison (Percentage)');
legend('Location', 'best', 'FontSize', 40)
xlim([data{1}.t(1), data{1}.t(end)]);
saveas(gcf, fullfile(figs_folder, 'P3 - Comparison Normalized Errors Percentage.png'));

fprintf('\nMax Normalized Error (%%):\n');
for i = 1:num_orders
    fprintf('Order n=%d: %.2f%%\n', data{i}.order, max_errors(i));
end
%% Comparison Accumulated Variance
figure('Units', 'normalized');
theme("light");
set(gcf, 'DefaultAxesFontSize', 24);
set(gcf, 'DefaultTextFontSize', 24);
hold on;
grid on;
for i = 1:num_orders
    plot(data{i}.t, data{i}.accumulated_variance, '-', 'LineWidth', 3, 'Color', data{i}.color, 'DisplayName', sprintf('n=%d', data{i}.order));
end
xlabel('Time [s]');
ylabel('Accumulated Variance [°C^2]');
title('Accumulated Variance: All Orders Overlaid');
legend('Location', 'best', 'FontSize', 40)
xlim([data{1}.t(1), data{1}.t(end)]);
saveas(gcf, fullfile(figs_folder, 'P3 - Comparison Accumulated Variance Overlaid.png'));

fprintf('Comparison plots saved to: %s\n', figs_folder);