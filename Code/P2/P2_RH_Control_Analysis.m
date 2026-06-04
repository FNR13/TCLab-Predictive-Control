%% =========================================================
%% P2 - Receding Horizon Control Analysis ECPD
%%
%% Group 17
%% Miguel Brás Simões - 99162
%% Henrique Miranda Póvoa - -102109
%% David de Jesus Alpalhão Marques - 100145
%% Ricardo Nascimento Francisco - 103093
%% =========================================================
clear; close all; clc;

%% Parameters
Q = 1;
B = 1;
C = 1;

% Representative R values
R_values = [0.1, 1, 10, 100];
R_labels  = {'R = 0.1', 'R = 1', 'R = 10', 'R = 100'};

% Color palette (distinct, colorblind-friendly)
colors = [
    0.0000, 0.4470, 0.7410;   % blue
    0.8500, 0.3250, 0.0980;   % orange
    0.4660, 0.6740, 0.1880;   % green
    0.7490, 0.1840, 0.5560;   % purple
];


Hmax     = 30;
H_values = 1:Hmax;


A_values    = [0.8, 1.2];
plant_names = {'A = 0.8  (planta estável em malha aberta)', ...
               'A = 1.2  (planta instável em malha aberta)'};
plant_short = {'A = 0.8', 'A = 1.2'};

%% =========================================================
%% FIGURE 1 & 2 — Open-loop 
%% =========================================================

x0      = 1;          % initial condition
k_steps = 0:20;

fig_ol = figure('Name', 'Análise em Malha Aberta', ...
                'Position', [50 50 1100 520]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:2
    A = A_values(p);
    y_free = x0 * A.^k_steps;

    nexttile;
    stem(k_steps, y_free, 'filled', 'LineWidth', 1.6, ...
        'Color', colors(p,:), 'MarkerFaceColor', colors(p,:));
    hold on;
    yline(0, 'k--', 'LineWidth', 1);
    yline(1, 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 1);

    grid on; box on;
    xlabel('Instante k', 'FontSize', 11);
    ylabel('y(k) = A^k \cdot x_0', 'FontSize', 11);
    title(sprintf('Resposta Livre em Malha Aberta — %s', plant_short{p}), ...
          'FontSize', 12, 'FontWeight', 'bold');

    if A < 1
        txt = sprintf('|A| = %.1f < 1\nSistema ESTÁVEL', A);
        stability_color = [0.1 0.6 0.1];
    else
        txt = sprintf('|A| = %.1f > 1\nSistema INSTÁVEL', A);
        stability_color = [0.8 0.1 0.1];
    end
    xlims = xlim; ylims = ylim;
    text(xlims(2)*0.97, ylims(2)*0.93, txt, ...
         'HorizontalAlignment', 'right', 'FontSize', 10.5, ...
         'Color', stability_color, 'FontWeight', 'bold', ...
         'BackgroundColor', [1 1 1 0.7], 'EdgeColor', stability_color, ...
         'Margin', 4);
    ylim([min(ylims(1), -0.1), max(ylims(2)*1.1, 1.15)]);
end

sgtitle('Análise de Estabilidade em Malha Aberta', ...
        'FontSize', 14, 'FontWeight', 'bold');


K_RH_all  = zeros(length(A_values), length(R_values), Hmax);
eig_RH_all = zeros(length(A_values), length(R_values), Hmax);
K_LQ_all  = zeros(length(A_values), length(R_values));
eig_LQ_all = zeros(length(A_values), length(R_values));

for p = 1:length(A_values)
    A = A_values(p);
    for r_idx = 1:length(R_values)
        R = R_values(r_idx);
        [K_LQ, ~, ~] = dlqr(A, B, Q, R);
        K_LQ_all(p, r_idx)  = K_LQ;
        eig_LQ_all(p, r_idx) = abs(A - B*K_LQ);

        for idx = 1:length(H_values)
            H = H_values(idx);
            [K_RH_all(p, r_idx, idx), ~] = receding_horizon_gain(A, B, C, R, H);
            eig_RH_all(p, r_idx, idx)    = abs(A - B*K_RH_all(p, r_idx, idx));
        end
    end
end

%% =========================================================
%% FIGURES 3 & 4 
%% =========================================================
for p = 1:length(A_values)
    A = A_values(p);

    fig = figure('Name', sprintf('Controlo RH — %s', plant_short{p}), ...
                 'Position', [100+p*30, 80, 1150, 680]);
    t   = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    
    nexttile;
    hold on;
    h_lines = gobjects(length(R_values), 1);
    for r_idx = 1:length(R_values)
        K_RH  = squeeze(K_RH_all(p, r_idx, :))';
        K_LQ  = K_LQ_all(p, r_idx);
        h_lines(r_idx) = plot(H_values, K_RH, 'o-', ...
            'Color', colors(r_idx,:), 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'MarkerFaceColor', colors(r_idx,:));
        yline(K_LQ, '--', 'Color', colors(r_idx,:), ...
              'LineWidth', 1.2, 'Alpha', 0.7);
    end
    grid on; box on;
    xlabel('Horizonte H', 'FontSize', 11);
    ylabel('Ganho K', 'FontSize', 11);
    title(sprintf('Convergência do Ganho K — Planta %s', plant_short{p}), ...
          'FontSize', 12, 'FontWeight', 'bold');

    leg_labels = cellfun(@(l) sprintf('K_{RH} (%s)', l), R_labels, ...
                         'UniformOutput', false);
    legend(h_lines, leg_labels, 'Location', 'best', 'FontSize', 9.5, ...
           'NumColumns', 2);
    xlim([1 Hmax]);

    annotation_str = '-- K_{LQ} (linha tracejada)';
    text(Hmax*0.98, min(ylim)*1.02 + range(ylim)*0.04, annotation_str, ...
         'HorizontalAlignment', 'right', 'FontSize', 8.5, ...
         'Color', [0.4 0.4 0.4]);

   
    nexttile;
    hold on;
    h_eig = gobjects(length(R_values), 1);
    for r_idx = 1:length(R_values)
        eig_RH = squeeze(eig_RH_all(p, r_idx, :))';
        h_eig(r_idx) = plot(H_values, eig_RH, 's-', ...
            'Color', colors(r_idx,:), 'LineWidth', 1.8, ...
            'MarkerSize', 5, 'MarkerFaceColor', colors(r_idx,:));
        yline(eig_LQ_all(p, r_idx), '--', 'Color', colors(r_idx,:), ...
              'LineWidth', 1.0, 'Alpha', 0.65);
    end
    hl = yline(1, 'k--', 'LineWidth', 2);
    text(Hmax*0.98, 1.02, 'Limite de Estabilidade  (|\lambda| = 1)', ...
         'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', [0.2 0.2 0.2]);

    grid on; box on;
    xlabel('Horizonte H', 'FontSize', 11);
    ylabel('|\lambda_{cl}|  =  |A - BK|', 'FontSize', 11);
    title(sprintf('Estabilidade: Magnitude do Valor Próprio |\lambda_{cl}| — Planta %s', plant_short{p}), ...
          'FontSize', 12, 'FontWeight', 'bold');

    leg_eig_labels = cellfun(@(l) sprintf('|\lambda| (%s)', l), R_labels, ...
                              'UniformOutput', false);
    legend(h_eig, leg_eig_labels, 'Location', 'best', 'FontSize', 9.5, ...
           'NumColumns', 2);
    xlim([1 Hmax]);

  
    yl = ylim;
    if yl(2) > 1
        fill([1 Hmax Hmax 1], [1 1 yl(2) yl(2)], [1 0.8 0.8], ...
             'FaceAlpha', 0.18, 'EdgeColor', 'none');
        text(Hmax*0.5, (1 + yl(2))/2, 'Região Instável', ...
             'HorizontalAlignment', 'center', 'FontSize', 9, ...
             'Color', [0.7 0.1 0.1], 'FontAngle', 'italic');
    end

    sgtitle(sprintf('Análise de Controlo em Horizonte Receding — %s', plant_short{p}), ...
            'FontSize', 13, 'FontWeight', 'bold');
end

%% =========================================================
%% FIGURE 5 — Summary: Minimum H to stabilise (A=1.2)
%% =========================================================
A = 1.2;
p = 2;   % index for A=1.2 in A_values

min_H_stable = NaN(1, length(R_values));
for r_idx = 1:length(R_values)
    eig_RH = squeeze(eig_RH_all(p, r_idx, :))';
    idx_stable = find(eig_RH < 1, 1, 'first');
    if ~isempty(idx_stable)
        min_H_stable(r_idx) = H_values(idx_stable);
    end
end

fig5 = figure('Name', 'Horizonte Mínimo de Estabilização — A=1.2', ...
              'Position', [200, 150, 700, 400]);

b = bar(1:length(R_values), min_H_stable, 0.55);
b.FaceColor = 'flat';
for r_idx = 1:length(R_values)
    b.CData(r_idx,:) = colors(r_idx,:);
    if ~isnan(min_H_stable(r_idx))
        text(r_idx, min_H_stable(r_idx) + 0.4, ...
             sprintf('H = %d', min_H_stable(r_idx)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10.5, ...
             'FontWeight', 'bold');
    else
        text(r_idx, 1, 'Não estabiliza', ...
             'HorizontalAlignment', 'center', 'FontSize', 9, ...
             'Color', 'r');
    end
end

grid on; box on;
xticks(1:length(R_values));
xticklabels(R_labels);
xlabel('Penalização de controlo R', 'FontSize', 11);
ylabel('Horizonte mínimo H para |\lambda_{cl}| < 1', 'FontSize', 11);
title({'Horizonte Mínimo de Estabilização — Planta A = 1.2', ...
       '(Planta instável em malha aberta)'}, ...
      'FontSize', 12, 'FontWeight', 'bold');
ylim([0, max(min_H_stable, [], 'omitnan') + 3]);

%% =========================================================
%% Console output
%% =========================================================
fprintf('\n=== Resumo dos Ganhos e Valores Próprios ===\n');
for p = 1:length(A_values)
    A = A_values(p);
    fprintf('\n--- %s ---\n', plant_names{p});
    fprintf('  %-12s %-12s %-15s %-18s %-15s\n', ...
            'R', 'K_LQ', '|lambda_LQ|', 'K_RH(H=30)', '|lambda_RH(H=30)|');
    for r_idx = 1:length(R_values)
        R    = R_values(r_idx);
        K_LQ = K_LQ_all(p, r_idx);
        eig_LQ = eig_LQ_all(p, r_idx);
        K_RH_final  = K_RH_all(p, r_idx, end);
        eig_RH_final = eig_RH_all(p, r_idx, end);
        fprintf('  %-12.3f %-12.6f %-15.6f %-18.6f %-15.6f\n', ...
                R, K_LQ, eig_LQ, K_RH_final, eig_RH_final);
    end
end

%% =========================================================
%% Local function: finite-horizon receding-horizon gain
%% =========================================================
function [K_RH, U_gain] = receding_horizon_gain(A, B, C, R, H)
    W  = zeros(H, H);
    Pi = zeros(H, 1);
    for row = 1:H
        Pi(row) = C * A^row;
        for col = 1:row
            W(row, col) = C * A^(row-col) * B;
        end
    end
    M      = W'*W + R*eye(H);
    U_gain = M \ (W'*Pi);
    K_RH   = U_gain(1);
end
