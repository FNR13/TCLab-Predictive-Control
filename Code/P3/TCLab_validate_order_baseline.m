% Validation comparison between low-order and selected TCLab models.
%
% This script uses the already saved state-space matrices in matfiles/model.
% It does not run a new system identification. The goal is to generate the
% extra validation figures requested for comparing n = 2 with the selected
% higher-order model.
%__________________________________________________________________________

clear
clc
close all

script_dir = fileparts(mfilename('fullpath'));
data_folder = fullfile(script_dir,'matfiles');
model_folder = fullfile(data_folder,'model');

baseline_order = 2;
selected_orders = [7 9];
orders_to_validate = unique([baseline_order selected_orders]);

validation_data = fullfile(data_folder,'openloop_data_2.mat');
mse_by_order = loadMSEByOrder(model_folder);

validation_results = struct([]);
for i = 1:numel(orders_to_validate)
    order = orders_to_validate(i);
    validation_results(i).order = order;
    validation_results(i).reference_mse = mse_by_order(order);
    validation_results(i).simulated_mse = plotValidationOrder(order, validation_data, model_folder, data_folder);
end

for i = 1:numel(selected_orders)
    plotValidationComparison(baseline_order, selected_orders(i), validation_data, model_folder, data_folder);
end

fprintf('Validation reference MSE values from MSE.mat:\n')
for i = 1:numel(validation_results)
    fprintf('n = %d: reference MSE = %.4f, simulated saved-model MSE = %.4f\n', ...
        validation_results(i).order, ...
        validation_results(i).reference_mse, ...
        validation_results(i).simulated_mse);
end

function mse_by_order = loadMSEByOrder(model_folder)
    mse_file = fullfile(model_folder,'MSE.mat');
    if exist(mse_file,'file')
        data = load(mse_file,'MSE');
        mse_by_order = data.MSE;
    else
        mse_by_order = nan(1,14);
    end
end

function mse = plotValidationOrder(order, validation_data, model_folder, output_folder)
    result = simulateValidation(order, validation_data, model_folder);
    mse = result.mse;

    fig = figure('Visible','off','Units','pixels','Position',[100 100 1200 750]);
    set(fig, 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', 16);
    set(fig, 'DefaultTextFontSize', 16);

    subplot(2,1,1), hold on, grid on
    plot(result.t,result.Dy,'.','MarkerSize',12,'Color',[0.10 0.35 0.70])
    plot(result.t,result.Dy_sim,'r--','LineWidth',2.5)
    title(sprintf('Model performance (n=%d) on validation dataset',order))
    xlabel('Time [s]')
    ylabel('\Delta y [deg C]')
    xlim([result.t(1),result.t(end)])
    lgd = legend('Experimental data','Model','Location','best');
    styleAxes()
    styleLegend(lgd)

    subplot(2,1,2), hold on, grid on
    plot(result.t,result.Dy-result.Dy_sim,'LineWidth',2.5,'Color',[0.15 0.15 0.15])
    xlabel('Time [s]')
    ylabel('\Delta y error [deg C]')
    xlim([result.t(1),result.t(end)])
    styleAxes()

    file_name = sprintf('P3 - SysId n%d - validation baseline.png',order);
    exportgraphics(fig,fullfile(output_folder,file_name),'Resolution',300)
    close(fig)
end

function plotValidationComparison(order_a, order_b, validation_data, model_folder, output_folder)
    result_a = simulateValidation(order_a, validation_data, model_folder);
    result_b = simulateValidation(order_b, validation_data, model_folder);

    fig = figure('Visible','off','Units','pixels','Position',[100 100 1200 750]);
    set(fig, 'Color', 'w');
    set(fig, 'DefaultAxesFontSize', 16);
    set(fig, 'DefaultTextFontSize', 16);

    subplot(2,1,1), hold on, grid on
    plot(result_a.t,result_a.Dy,'.','MarkerSize',12,'Color',[0.10 0.35 0.70])
    plot(result_a.t,result_a.Dy_sim,'--','LineWidth',2.5,'Color',[0.85 0.20 0.20])
    plot(result_b.t,result_b.Dy_sim,'-','LineWidth',2.5,'Color',[0.10 0.55 0.25])
    title(sprintf('Validation comparison: n=%d vs n=%d',order_a,order_b))
    xlabel('Time [s]')
    ylabel('\Delta y [deg C]')
    xlim([result_a.t(1),result_a.t(end)])
    lgd = legend('Experimental data',sprintf('Model n=%d',order_a),sprintf('Model n=%d',order_b),'Location','best');
    styleAxes()
    styleLegend(lgd)

    subplot(2,1,2), hold on, grid on
    plot(result_a.t,result_a.Dy-result_a.Dy_sim,'--','LineWidth',2.5,'Color',[0.85 0.20 0.20])
    plot(result_b.t,result_b.Dy-result_b.Dy_sim,'-','LineWidth',2.5,'Color',[0.10 0.55 0.25])
    xlabel('Time [s]')
    ylabel('\Delta y error [deg C]')
    xlim([result_a.t(1),result_a.t(end)])
    lgd = legend(sprintf('Error n=%d',order_a),sprintf('Error n=%d',order_b),'Location','best');
    styleAxes()
    styleLegend(lgd)

    file_name = sprintf('P3 - SysId n%d vs n%d - validation.png',order_a,order_b);
    exportgraphics(fig,fullfile(output_folder,file_name),'Resolution',300)
    close(fig)
end

function result = simulateValidation(order, validation_data, model_folder)
    model_file = fullfile(model_folder,sprintf('singleheater_model_%d.mat',order));
    model = load(model_file,'A','B','C','y_ss','u_ss');
    data = load(validation_data,'y','u','t');

    t = double(data.t(:))';
    y = double(data.y(1,:));
    u = double(data.u(1,:));
    Dy = y - model.y_ss;
    Du = u - model.u_ss;

    Dy_sim = simulateSavedModel(model.A, model.B, model.C, Du, Dy);

    result.order = order;
    result.t = t;
    result.Dy = Dy;
    result.Du = Du;
    result.Dy_sim = Dy_sim;
    result.mse = mean((Dy_sim - Dy).^2);
end

function Dy_sim = simulateSavedModel(A, B, C, Du, Dy)
    n = size(A,1);
    N = numel(Du);

    observability_rows = zeros(N-1,n);
    known_response = zeros(N-1,1);

    A_power = eye(n);
    forced_state = zeros(n,1);
    for k = 2:N
        forced_state = A*forced_state + B*Du(k-1);
        A_power = A*A_power;
        observability_rows(k-1,:) = C*A_power;
        known_response(k-1) = C*forced_state;
    end

    initial_state = observability_rows \ (Dy(2:end)' - known_response);

    state = nan(n,N);
    Dy_sim = nan(1,N);
    state(:,1) = initial_state;
    Dy_sim(1) = Dy(1);

    for k = 1:N-1
        state(:,k+1) = A*state(:,k) + B*Du(k);
        Dy_sim(k+1) = C*state(:,k+1);
    end
end

function styleAxes()
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.82 0.82 0.82];
    ax.GridAlpha = 0.7;
    ax.Title.Color = 'k';
    ax.XLabel.Color = 'k';
    ax.YLabel.Color = 'k';
end

function styleLegend(lgd)
    lgd.Color = 'w';
    lgd.TextColor = 'k';
    lgd.EdgeColor = [0.45 0.45 0.45];
end
