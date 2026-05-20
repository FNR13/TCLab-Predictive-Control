% P2 - Receding-horizon (RH) vs infinite-horizon LQ gain.
% Reproduces figures 3, 6 and 7 of the project report (Portuguese labels).
%
% Plants:  x(k+1) = A x(k) + u(k),    B = 1, C = 1, Q = 1.
%   A = 1.2  (open-loop unstable),  A = 0.8 (open-loop stable).
% Sweep:   R in {0.1, 1, 10, 100},  H = 1..30.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
out_dir    = fullfile(script_dir, '..', '..', '..', 'Imagens');
if ~exist(out_dir,'dir'); mkdir(out_dir); end

set(0,'DefaultFigureColor','w','DefaultAxesColor','w', ...
      'DefaultAxesXColor','k','DefaultAxesYColor','k', ...
      'DefaultAxesGridColor',[0.15 0.15 0.15],'DefaultTextColor','k', ...
      'DefaultAxesFontSize',11,'DefaultLegendColor','w', ...
      'DefaultLegendTextColor','k','DefaultLegendEdgeColor','k', ...
      'DefaultFigureInvertHardcopy','off');
save_white = @(fig, name) exportgraphics(fig, fullfile(out_dir,name), ...
                                         'BackgroundColor','white','Resolution',150);

B  = 1;  C = 1;  Q = 1;
R_list = [0.1, 1, 10, 100];
H_MAX  = 30;
H_grid = 1:H_MAX;

% Colour palette matching the report
%   R = 0.1   -> blue
%   R = 1     -> orange
%   R = 10    -> green
%   R = 100   -> magenta
colors = [ ...
    0.000, 0.447, 0.741;
    0.850, 0.325, 0.098;
    0.466, 0.674, 0.188;
    0.929, 0.000, 0.541];

% ========================================================================
% Figure 3 - Open-loop stability
% ========================================================================
K_free = 21;
k = 0:K_free-1;
y_stable   = (0.8).^k;
y_unstable = (1.2).^k;

fig = figure('Position',[100 100 1050 450],'Color','w');
sgtitle('Análise de Estabilidade em Malha Aberta','FontWeight','bold','Color','k');

subplot(1,2,1);
stem(k, y_stable, 'filled','LineWidth',1.0, ...
     'Color',colors(1,:),'MarkerFaceColor',colors(1,:));
grid on; box on;
title('Resposta Livre em Malha Aberta — A = 0.8');
xlabel('Instante k'); ylabel('y(k) = A^k x_0');
ylim([0, 1.05]); xlim([0, 20]);
text(11, 0.92, sprintf('|A| = 0.8 < 1\nSistema ESTÁVEL'), ...
     'Color',[0 0.5 0],'FontWeight','bold','HorizontalAlignment','left', ...
     'BackgroundColor',[0.85 1 0.85],'EdgeColor',[0 0.5 0]);

subplot(1,2,2);
stem(k, y_unstable, 'filled','LineWidth',1.0, ...
     'Color',colors(2,:),'MarkerFaceColor',colors(2,:));
grid on; box on;
title('Resposta Livre em Malha Aberta — A = 1.2');
xlabel('Instante k'); ylabel('y(k) = A^k x_0');
xlim([0, 20]); ylim([0, 1.05*max(y_unstable)]);
text(2, max(y_unstable)*0.92, sprintf('|A| = 1.2 > 1\nSistema INSTÁVEL'), ...
     'Color',[0.7 0 0],'FontWeight','bold','HorizontalAlignment','left', ...
     'BackgroundColor',[1 0.88 0.88],'EdgeColor',[0.7 0 0]);

save_white(fig,'P2_openloop_stability.png'); close(fig);
fprintf('P2_openloop_stability.png written\n');

% ========================================================================
% Sweep + plots for each plant
% ========================================================================
A_list = [1.2, 0.8];
A_tags = {'1.2','0.8'};
fname  = {'P2_RH_unstable.png','P2_RH_stable.png'};

fprintf('\n--- P2: receding-horizon vs LQ ---\n');
fprintf('%-8s %-8s %-12s %-12s %-14s\n', ...
        'A','R','K_LQ','K_RH(H_max)','|A-B*K_RH|');

for ia = 1:numel(A_list)
    A = A_list(ia);

    K_LQ_vec = zeros(1, numel(R_list));
    lam_LQ   = zeros(1, numel(R_list));
    K_RH_mat = zeros(numel(R_list), H_MAX);
    lam_RH   = zeros(numel(R_list), H_MAX);

    for ir = 1:numel(R_list)
        R = R_list(ir);
        K_LQ = dlqr(A, B, Q, R);
        K_LQ_vec(ir) = K_LQ;
        lam_LQ(ir)   = abs(A - B*K_LQ);

        for ih = 1:H_MAX
            H  = ih;
            W  = zeros(H);  Pi = zeros(H,1);
            for i = 1:H
                Pi(i) = C * A^i;
                for j = 1:i
                    W(i,j) = C * A^(i-j) * B;
                end
            end
            M    = W.'*W + R*eye(H);
            K_RH = [1, zeros(1,H-1)] * (M \ (W.' * Pi));
            K_RH_mat(ir,ih) = K_RH;
            lam_RH(ir,ih)   = abs(A - B*K_RH);
        end
        fprintf('%-8.2f %-8g %-12.4f %-12.4f %-14.4f\n', ...
                A, R, K_LQ, K_RH_mat(ir,end), lam_RH(ir,end));
    end

    % --- figure ----------------------------------------------------------
    fig = figure('Position',[100 100 950 720],'Color','w');
    sgtitle(sprintf('Análise de Controlo em Horizonte Receding — A = %s', A_tags{ia}), ...
            'FontWeight','bold','Color','k');

    % ---- top panel: gains
    subplot(2,1,1); hold on; grid on; box on;
    title(sprintf('Convergência do Ganho K — Planta A = %s', A_tags{ia}));
    for ir = 1:numel(R_list)
        plot(H_grid, K_RH_mat(ir,:), 'o-', 'Color', colors(ir,:), ...
             'MarkerSize',4, 'MarkerFaceColor', colors(ir,:), 'LineWidth', 1.2, ...
             'DisplayName', sprintf('K_{RH} (R = %g)', R_list(ir)));
        yline(K_LQ_vec(ir), '--', 'Color', colors(ir,:), 'LineWidth', 0.9, ...
              'HandleVisibility','off');
    end
    xlabel('Horizonte H'); ylabel('Ganho K');
    legend('Location','best','NumColumns',2);
    xlim([1, H_MAX]);
    ylim_top = ylim;
    ylim([max(0, ylim_top(1)), ylim_top(2)*1.05]);

    % ---- bottom panel: closed-loop eigenvalue magnitude
    subplot(2,1,2); hold on; grid on; box on;
    title('Estabilidade: Magnitude do Valor Próprio |\lambda|');

    lam_max = max(lam_RH(:));
    if lam_max >= 1
        ymax = max(1.25, lam_max*1.05);
        % shade unstable region
        patch([H_grid(1) H_grid(end) H_grid(end) H_grid(1)], ...
              [1 1 ymax ymax], [1 0.85 0.85], ...
              'EdgeColor','none','FaceAlpha',0.55, ...
              'HandleVisibility','off');
        text(H_MAX*0.5, (1+ymax)/2, 'Região Instável', ...
             'Color',[0.6 0 0],'FontWeight','bold','HorizontalAlignment','center', ...
             'FontAngle','italic');
    else
        ymax = 1.05;
    end

    for ir = 1:numel(R_list)
        plot(H_grid, lam_RH(ir,:), 'o-', 'Color', colors(ir,:), ...
             'MarkerSize',4, 'MarkerFaceColor', colors(ir,:), 'LineWidth', 1.2, ...
             'DisplayName', sprintf('R = %g', R_list(ir)));
    end
    yline(1, 'k--', 'Limite de Estabilidade |\lambda| = 1', ...
          'LineWidth', 1.0, 'LabelHorizontalAlignment','right', ...
          'HandleVisibility','off');
    xlabel('Horizonte H'); ylabel('|\lambda_{cl}| = |A - B K|');
    legend('Location','best','NumColumns',2);
    xlim([1, H_MAX]); ylim([0, ymax]);

    save_white(fig, fname{ia}); close(fig);
    fprintf('%s written\n', fname{ia});
end
