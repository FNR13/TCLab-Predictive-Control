% Identification of TCLab model for single heater
%
% Loads dataset 'openloop_data_1.mat' and identifies a discrete linear
% time-invariant model for the incremental dynamics around an equilibrium
% point. Validates the model for a different dataset 'openloop_data_2.mat'.
%
% Functions called: ssest, findstates.
%
% Afonso Botelho and J. Miranda Lemos, IST, May 2023
%__________________________________________________________________________

%% Initialization
clear
clc
close all

number_of_iterations = 14;
selected_model_order = 9;
selected_models_to_save_image = [1, 2, 6, 9, 14];

MSE = zeros(1, number_of_iterations);
Variance = zeros(1, number_of_iterations);

script_dir = fileparts(mfilename('fullpath'));
common_dir = fullfile(script_dir, '..', 'Common');
figs_folder = fullfile(script_dir, '..', 'figs');
data_folder = fullfile(script_dir,'matfiles');
model_folder = fullfile(data_folder,'model');

%% Identification
for n=1:number_of_iterations
    % Load data and select the output/input for the first heater only
    load(fullfile(data_folder,'openloop_data_1.mat'),'y','u','t');
    u = u(1,:);
    y = y(1,:);
    
    % Choose interval for initial equilibrium
    k_ss_begin = 201; % initial sample
    k_ss_end = 400; % final sample
    
    % Compute steady-state output/input from initial equilibrium
    y_ss = mean(y(:,k_ss_begin:k_ss_end),2);
    u_ss = u(:,k_ss_begin);
    
    % Truncate initial transient
    t = t(k_ss_begin:end-1);
    u = u(:,k_ss_begin:end-1);
    y = y(:,k_ss_begin:end-1);
    
    % Compute incremental output/input
    Dy = y - y_ss;
    Du = u - u_ss;
    
    % Identify state-space system for incremental dynamics
    % n = 1;
    Ts = t(2) - t(1);
    sys = ssest(Du',Dy',n,'Ts',Ts);
    [A,B,C,~,Ke] = idssdata(sys);
    e_var = sys.NoiseVariance;

    save(fullfile(model_folder,sprintf('singleheater_model_%d.mat',n)),'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
    if n == selected_model_order
        save(fullfile(model_folder,'singleheater_model.mat'),'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
        % Also publish the selected model to Common/ so that the P4 and
        % P5 scripts find it on a fixed path.
        save(fullfile(common_dir,'singleheater_model.mat'),'A','B','C','Ke','e_var','y_ss','u_ss','Ts');
    end
    
    %% Test on dataset 1, with which the model was identified
    
    % Initializations
    N = length(t);
    Dy_sim = nan(1,N);
    Dx_sim = nan(n,N);
    
    % Find initial incremental state that best fits the data given the identified model
    Dx0 = findstates(sys,iddata(Dy',Du',Ts));
    
    % Set initial conditions
    Dy_sim(:,1) = Dy(:,1);
    Dx_sim(:,1) = Dx0;
    
    % Propagate model
    for k = 1:N-1
        Dx_sim(:,k+1) = A*Dx_sim(:,k) + B*Du(:,k);
        Dy_sim(:,k+1) = C*Dx_sim(:,k+1);
    end
    
    % Plot results
    figure('Units','normalized');
    theme("light")
    set(gcf, 'DefaultAxesFontSize', 24) 
    set(gcf, 'DefaultTextFontSize', 24)
    subplot(2,1,1), hold on, grid on   
    title(sprintf('Model performance (n=%d) on identification dataset',n))
    width = 5;
    plot(t,Dy,'.','MarkerSize',15)
    plot(t,Dy_sim,'r--', 'LineWidth',width-2)
    xlabel('Time [s]')
    ylabel('\Delta y [°C]')
    xlim([t(1),t(end)]);
    legend('Experimental data','Model','Location','best');
    subplot(2,1,2), hold on, grid on   
    stairs(t,Du,'LineWidth',width)
    xlabel('Time [s]')
    ylabel('\Delta u [%]')
    xlim([t(1),t(end)]);
    
    if ismember(n, selected_models_to_save_image)
        filename_training = sprintf('P3 - SysId n%d - training.png', n);
        saveas(gcf, fullfile(figs_folder, filename_training));
    end
    %% Test on dataset 2, with which the model was not identified

    % Load data and select the output/input for the first heater only
    load(fullfile(data_folder,'openloop_data_2.mat'),'y','u','t');
    u = u(1,:);
    y = y(1,:);
    
    % Compute incremental output/input
    Dy2 = y - y_ss;
    Du2 = u - u_ss;
    
    % Initializations
    N = length(t);
    Dy2_sim = nan(1,N);
    Dx2_sim = nan(n,N);
    
    % Find initial incremental state that best fits the data given the identified model
    Dx02 = findstates(sys,iddata(Dy2',Du2',Ts));
    
    % Set initial conditions
    Dy2_sim(:,1) = Dy2(:,1);
    Dx2_sim(:,1) = Dx02;
    
    % Propagate model
    for k = 1:N-1
        Dx2_sim(:,k+1) = A*Dx2_sim(:,k) + B*Du2(:,k);
        Dy2_sim(:,k+1) = C*Dx2_sim(:,k+1);
    end
    
    % Plot results
    figure('Units','normalized')
    theme("light")
    set(gcf, 'DefaultAxesFontSize', 24) 
    set(gcf, 'DefaultTextFontSize', 24)
    subplot(2,1,1), hold on, grid on   
    title(sprintf('Model performance (n=%d) on validation dataset',n))
    plot(t,Dy2,'.','MarkerSize',15)
    plot(t,Dy2_sim,'r--', 'LineWidth',width-2)
    xlabel('Time [s]')
    ylabel('\Delta{y} [°C]')
    xlim([t(1),t(end)]);
    legend('Experimental data','Model','Location','best');
    subplot(2,1,2), hold on, grid on   
    plot(t,Dy2-Dy2_sim,'MarkerSize',5, 'LineWidth', width)
    xlabel('Time [s]')
    ylabel('\Delta{y} Error [°C]')
    xlim([t(1),t(end)]);
    
    if ismember(n, selected_models_to_save_image)
        filename_validation = sprintf('P3 - SysId n%d - validation.png', n);
        saveas(gcf, fullfile(figs_folder, filename_validation));
    end
    
    error = Dy2_sim-Dy2;
    mse = sum(error.^2)/N;
    variance = var(error);
    
    accumulated_variance = zeros(1, N);
    for k = 1:N
        accumulated_variance(k) = var(error(1:k)); % Variance up to time k
    end

    if ismember(n, selected_models_to_save_image)
        figure('Units','normalized');
        theme("light")
        set(gcf, 'DefaultAxesFontSize', 24);
        set(gcf, 'DefaultTextFontSize', 24);
        plot(t, accumulated_variance, 'b-', 'LineWidth', 2);
        xlabel('Time [s]');
        ylabel('Accumulated Variance [°C^2]');
        title(sprintf('Accumulated Error Variance (n=%d)', n));
        grid on;
        saveas(gcf, fullfile(figs_folder, sprintf('P3 - SysId n%d - accumulated_variance.png', n)));
    end

    MSE(n) = mse;
    Variance(n) = variance;
    fprintf('Error between propagated and measured output for %.0f state on: MSE: %.4f Variance %4f\n',n, mse, variance);

    save(...
        fullfile(model_folder, sprintf('simulation_results_n%d.mat', n)), ...
        't', 'Dy2', 'Dy2_sim', 'error', 'accumulated_variance', 'Du2', 'Dy', 'Dy_sim', 'Du', 'y_ss', 'u_ss', 'Ts' ...
    );
end

%% Erro Plotting
% Create MSE bar graph
figure;
theme("light")
set(gcf, 'DefaultAxesFontSize', 24) 
set(gcf, 'DefaultTextFontSize', 24)
bar(MSE, 'FaceColor', [0.529, 0.808, 0.980]); % Sky blue color
xlabel('Order n');
ylabel('MSE');
title('MSE by Order n');
grid on;

hold on;
for i = 1:length(MSE)
    % Get the x and y position for the text
    x_pos = i;
    y_pos = MSE(i);

    % Add text with the MSE value, slightly above the bar
    text(x_pos, y_pos, sprintf('%.4f', MSE(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 16, ...
         'Color', 'black');
end
hold off;

filename_MSE = 'P3 - SysId Error - MSE by Order n.png';
saveas(gcf, fullfile(figs_folder, filename_MSE));

% Create Variance bar graph
figure;
theme("light")
set(gcf, 'DefaultAxesFontSize', 24) 
set(gcf, 'DefaultTextFontSize', 24)
bar(Variance, 'FaceColor', [0.529, 0.808, 0.980]); % Sky blue color
xlabel('Order n');
ylabel('Variance');
title('Variance by Order n');
grid on;

hold on;
for i = 1:length(Variance)
    % Get the x and y position for the text
    x_pos = i;
    y_pos = Variance(i);

    % Add text with the MSE value, slightly above the bar
    text(x_pos, y_pos, sprintf('%.4f', Variance(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'bottom', ...
         'FontSize', 16, ...
         'Color', 'black');
end
hold off;

filename_Variance = 'P3 - SysId Error - Variance by Order n.png';
saveas(gcf, fullfile(figs_folder, filename_Variance));
fprintf('\n')
save(fullfile(model_folder,'MSE.mat'), 'MSE')
save(fullfile(model_folder,'Variance .mat'), 'Variance')
fprintf('Selected model n = %d saved as singleheater_model.mat\n',selected_model_order);
%--------------------------------------------------------------------------
% End of File
