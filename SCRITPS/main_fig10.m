clc; clear; close all;
rng(1);

%%
config;
configplot;

%% =========================================================
% 真值目标 PSM
vec = S(2, :);
k   = vec(1);
rs  = vec(2);
phx = vec(3);
phv = vec(4);

S = [1,                  k  * exp(1j * phx);
     k * exp(1j * phx),  rs * exp(1j * phv)];

%% 真值干扰极化
j_gamma = pi / 3;
j_eta   = pi / 2;
pj_true = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];
pj_true = pj_true / norm(pj_true);

%% =========================================================
% 参数设置
polar_num = 30;
train_len = N_temp;

% 这里只取两个 JSR，作为主文后续画 6 条线的候选
JSR_fix = [0, 10, 20];

% 注入的 pj 子空间误差：
% eps_p = 1 - |phat' * pj_true|^2 = sin^2(theta)
eps_p_list = [0, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2];
theta_deg  = asind(sqrt(eps_p_list));

MC_noise = 100;   % 噪声 Monte-Carlo
MC_eps   = 30;    % 每个 eps 下，随机相位扰动 Monte-Carlo

%% 5 个轨道都跑
track2s = [14, 2, 3, 4, 11];
name_map = containers.Map('KeyType','double','ValueType','char');
name_map(14) = 'LIN';
name_map(2)  = 'HV';
name_map(3)  = 'RCL';
name_map(4)  = 'GND';
name_map(11) = 'OUR';

names_sel = cell(1, numel(track2s));
for i = 1:numel(track2s)
    names_sel{i} = name_map(track2s(i));
end

%% =========================================================
% 结果数组
% size = [num_track, num_JSR, num_eps]
sinr_out = zeros(numel(track2s), numel(JSR_fix), numel(eps_p_list));
gain_out = zeros(numel(track2s), numel(JSR_fix), numel(eps_p_list));

%% =========================================================
% 主循环
for it = 1:numel(track2s)

    track2 = track2s(it);
    fprintf('TRACK %2d / %2d : %s\n', it, numel(track2s), names_sel{it});

    for ij = 1:numel(JSR_fix)

        jsr_now = JSR_fix(ij);
        fprintf('  JSR = %+3d dB\n', jsr_now);

        for ie = 1:numel(eps_p_list)

            eps_now = eps_p_list(ie);

            s_tmp = zeros(1, MC_eps);
            g_tmp = zeros(1, MC_eps);

            for im = 1:MC_eps

                % 直接构造带误差的 phat，不经过 pulse-1
                phat = inject_pj_subspace_error(pj_true, eps_now);

                [s_tmp(im), g_tmp(im)] = run_pulse2_track_once( ...
                    jsr_now, polar_num, track2, S, ...
                    train_len, pj_true, phat, MC_noise);
            end

            sinr_out(it, ij, ie) = mean(s_tmp);
            gain_out(it, ij, ie) = mean(g_tmp);

            fprintf('    eps_p = %.4g | theta = %.3f deg | SINR = %.3f dB | Gain = %.3f dB\n', ...
                eps_now, theta_deg(ie), sinr_out(it, ij, ie), gain_out(it, ij, ie));
        end
    end
end

%% =========================================================
% 保存数据
save('z_tsp1_fig10.mat', ...
    'sinr_out', 'gain_out', ...
    'eps_p_list', 'theta_deg', 'JSR_fix', ...
    'track2s', 'names_sel', ...
    'MC_noise', 'MC_eps');

fprintf('\nSaved to z_tsp1_fig10.mat\n');

%% =========================================================
% 预览图 1：每个 JSR 下，5 个轨道的 SINR 曲线
for ij = 1:numel(JSR_fix)

    figure('Color', 'w'); hold on; grid on; box on;

    for it = 1:numel(track2s)
        plot(theta_deg, squeeze(sinr_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{it}, ...
            'DisplayName', names_sel{it});
    end

    xlabel('Principal-angle error in $\hat{\mathbf p}_j$ (deg)', ...
        'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Output SINR (dB)', ...
        'Interpreter', 'latex', 'FontSize', 18);
    title(sprintf('Pulse-2 output SINR vs $\\hat{p}_j$ error (JSR = %d dB)', JSR_fix(ij)), ...
        'Interpreter', 'latex', 'FontSize', 16);

    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex');

    set(gcf, 'Position', [200 160 800 600]);
end

%% 预览图 2：每个 JSR 下，5 个轨道的 GAIN 曲线
for ij = 1:numel(JSR_fix)

    figure('Color', 'w'); hold on; grid on; box on;

    for it = 1:numel(track2s)
        plot(theta_deg, squeeze(gain_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{it}, ...
            'DisplayName', names_sel{it});
    end

    xlabel('Principal-angle error in $\hat{\mathbf p}_j$ (deg)', ...
        'Interpreter', 'latex', 'FontSize', 18);
    ylabel('SINR gain (dB)', ...
        'Interpreter', 'latex', 'FontSize', 18);
    title(sprintf('Pulse-2 SINR gain vs $\\hat{p}_j$ error (JSR = %d dB)', JSR_fix(ij)), ...
        'Interpreter', 'latex', 'FontSize', 16);

    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontName', 'Times New Roman', ...
        'FontSize', 16, ...
        'LineWidth', 1.0, ...
        'TickLabelInterpreter', 'latex');

    set(gcf, 'Position', [260 180 800 600]);
end

%% =========================================================
% 预览图 3：固定 5 个轨道，两个 JSR 叠在一张图里（SINR）
% 这个图就是你后面主文可能会从里面筛 3 个轨道 * 2 个 JSR = 6 条线 的母图
figure('Color', 'w'); hold on; grid on; box on;

cc = 1;
for it = 1:numel(track2s)
    for ij = 1:numel(JSR_fix)
        this_name = sprintf('%s, JSR=%d dB', names_sel{it}, JSR_fix(ij));
        plot(theta_deg, squeeze(sinr_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{ij}, ...
            'DisplayName', this_name);
        cc = cc + 1;
    end
end

xlabel('Principal-angle error in $\hat{\mathbf p}_j$ (deg)', ...
    'Interpreter', 'latex', 'FontSize', 18);
ylabel('Output SINR (dB)', ...
    'Interpreter', 'latex', 'FontSize', 18);
title('All candidate curves for Pulse-2 robustness (SINR)', ...
    'Interpreter', 'latex', 'FontSize', 16);

legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'LineWidth', 1.0, ...
    'TickLabelInterpreter', 'latex');

set(gcf, 'Position', [120 120 1000 620]);

%% 预览图 4：固定 5 个轨道，两个 JSR 叠在一张图里（GAIN）
figure('Color', 'w'); hold on; grid on; box on;

cc = 1;
for it = 1:numel(track2s)
    for ij = 1:numel(JSR_fix)
        this_name = sprintf('%s, JSR=%d dB', names_sel{it}, JSR_fix(ij));
        plot(theta_deg, squeeze(gain_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{ij}, ...
            'DisplayName', this_name);
        cc = cc + 1;
    end
end

xlabel('Principal-angle error in $\hat{\mathbf p}_j$ (deg)', ...
    'Interpreter', 'latex', 'FontSize', 18);
ylabel('SINR gain (dB)', ...
    'Interpreter', 'latex', 'FontSize', 18);
title('All candidate curves for Pulse-2 robustness (Gain)', ...
    'Interpreter', 'latex', 'FontSize', 16);

legend('Location', 'eastoutside', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'FontName', 'Times New Roman', ...
    'FontSize', 16, ...
    'LineWidth', 1.0, ...
    'TickLabelInterpreter', 'latex');

set(gcf, 'Position', [140 140 1000 620]);

%% =========================================================
function [sinr_avg, gain_avg] = run_pulse2_track_once( ...
    jsr, polar_num, track2, S, train_len, pj_true, phat, MC)

    config;

    SNR = 10;
    Ps  = 1 * 10^(SNR / 10);
    Pi  = 10^((jsr + SNR) / 10);

    lambda_scale = 5e-2;

    win       = start_idx : end_idx;
    train_idx = start_idx : (start_idx + train_len - 1);

    %% 只做 pulse-2
    opts = struct();
    opts.opt2_pj = phat;

    pk2 = gen_pk(polar_num, track2, S, opts);

    [lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);
    [lfm_h2, lfm_v2, ~] = tx_polarize(lfm_temp, pk2);

    [lfm_scatter2, ~, ~] = build_target_echo(S, lfm_h2, lfm_v2, N_PRI, delay_samples);
    t = sum(lfm_scatter2(:, win), 2);
    if norm(t) == 0
        t = [1; 0];
    end

    % 干扰始终用真实 pj_true 生成
    smsp2 = gen_smsp_polar(lfm_temp, pk2, pj_true, T, B, fs, fc, ...
        M_smsp, delta_smsp, N_PRI, delay_samples, Pi);

    s = zeros(1, MC);
    g = zeros(1, MC);

    for ii = 1:MC

        noise = [sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                 sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];

        rsig2 = lfm_scatter2 + smsp2 + noise;

        sinr_pre = get_sinr(rsig2, lfm_scatter2);

        Ytr2 = rsig2(:, train_idx);
        Rhat = (Ytr2 * Ytr2') / max(1, size(Ytr2, 2));
        Rhat = (Rhat + Rhat') / 2;

        lambda = lambda_scale * trace(Rhat) / 2;
        Rdl    = Rhat + lambda * eye(2);

        C = [t, phat];
        f = [1; 0];

        Gm = C' * (Rdl \ C);
        if rcond(Gm) < 1e-6
            Gm = Gm + 1e-6 * eye(2);
        end

        w_lcmv = (Rdl \ C) * (Gm \ f);

        s(ii) = get_sinr(w_lcmv' * rsig2, w_lcmv' * lfm_scatter2);
        g(ii) = s(ii) - sinr_pre;
    end

    sinr_avg = mean(s);
    gain_avg = mean(g);
end

%% =========================================================
function phat = inject_pj_subspace_error(pj_true, eps_p)
% 注入受控子空间误差：
% eps_p = 1 - |phat' * pj_true|^2 = sin^2(theta)

    pj_true = pj_true / max(norm(pj_true), 1e-12);

    if eps_p <= 0
        phat = pj_true;
        return;
    end

    theta = asin(sqrt(eps_p));

    p_perp = [-conj(pj_true(2)); conj(pj_true(1))];
    p_perp = p_perp / max(norm(p_perp), 1e-12);

    phi = 2 * pi * rand;

    phat = cos(theta) * pj_true + exp(1j * phi) * sin(theta) * p_perp;
    phat = phat / max(norm(phat), 1e-12);
end