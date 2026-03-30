clc; clear; close all;
rng(1);

%%
config;
configplot;

load psm.mat;

%% =========================================================
% True target PSM
vec = s1(2, :);
k   = vec(1);
rs  = vec(2);
phx = vec(3);
phv = vec(4);

S = [1,                  k  * exp(1j * phx);
     k * exp(1j * phx),  rs * exp(1j * phv)];

%% True jammer polarization
j_gamma = pi / 3;
j_eta   = pi / 2;
pj_true = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];
pj_true = pj_true / norm(pj_true);

%% =========================================================
% Parameters
polar_num = 10;                 % consistent with Table I / reviewer concern
train_len = N_temp;

% two representative JSRs
JSR_fix = [20];

% phase mismatch bound (deg)
phi_list = [0, 1, 2, 3, 5, 7, 10];

MC_noise = 100;                 % noise Monte-Carlo
MC_phase = 20;                  % outer Monte-Carlo over random phase errors

%% Five trajectories
track1s = [14, 2, 3, 4, 10];
track2s = [14, 2, 3, 4, 11];
names_sel = {'LIN', 'HV', 'RCL', 'GND', 'OUR'};

%% =========================================================
% result arrays
% size = [num_track, num_JSR, num_phi]
sinr_out = zeros(numel(track1s), numel(JSR_fix), numel(phi_list));
gain_out = zeros(numel(track1s), numel(JSR_fix), numel(phi_list));

%% =========================================================
% main loop
for it = 1:numel(track1s)

    fprintf('TRACK %2d / %2d : %s\n', it, numel(track1s), names_sel{it});

    track1 = track1s(it);
    track2 = track2s(it);

    for ij = 1:numel(JSR_fix)

        jsr_now = JSR_fix(ij);
        fprintf('  JSR = %+3d dB\n', jsr_now);

        for ip = 1:numel(phi_list)

            phi_now = phi_list(ip);

            s_tmp = zeros(1, MC_phase);
            g_tmp = zeros(1, MC_phase);

            for im = 1:MC_phase
                [s_tmp(im), g_tmp(im)] = run_once_phase_mismatch( ...
                    jsr_now, polar_num, track1, track2, S, ...
                    train_len, j_gamma, j_eta, ...
                    phi_now, MC_noise);
            end

            sinr_out(it, ij, ip) = mean(s_tmp);
            gain_out(it, ij, ip) = mean(g_tmp);

            fprintf('    phi = %2d deg | SINR = %.3f dB | Gain = %.3f dB\n', ...
                phi_now, sinr_out(it, ij, ip), gain_out(it, ij, ip));
        end
    end
end

%% =========================================================
% save
save('z_tsp1_fig12.mat', ...
    'sinr_out', 'gain_out', ...
    'phi_list', 'JSR_fix', ...
    'track1s', 'track2s', 'names_sel', ...
    'MC_noise', 'MC_phase');

fprintf('\nSaved to z_tsp1_fig12.mat\n');

%% =========================================================
% Preview: gain curves
for ij = 1:numel(JSR_fix)

    figure('position', p); hold on; grid on; box on;

    for it = 1:numel(track1s)
        plot(phi_list, squeeze(gain_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{it}, ...
            'DisplayName', names_sel{it});
    end

    xlabel('Phase mismatch bound (deg)', 'Interpreter', 'latex');
    ylabel('SINR Gain (dB)', 'Interpreter', 'latex');

    legend('Location', 'southeast', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Times New Roman', ...
        'TickLabelInterpreter', 'latex');
    ylim([12 20])

end

%% =========================================================
function [sinr_avg, gain_avg] = run_once_phase_mismatch( ...
    jsr, polar_num, track1, track2, S, ...
    train_len, j_gamma, j_eta, phi_deg, MC)

    config;

    SNR = 10;
    Ps  = 1 * 10^(SNR / 10);
    Pi  = 10^((jsr + SNR) / 10);

    lambda_scale = 5e-2;

    win       = start_idx : end_idx;
    train_idx = start_idx : (start_idx + train_len - 1);

    %% true jammer polarization
    pj = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];
    pj = pj / max(norm(pj), 1e-12);

    %% ---------- Pulse 1: ideal design, mismatched execution ----------
    pk1_ideal = gen_pk(polar_num, track1, S);
    pk1_exec  = apply_phase_mismatch(pk1_ideal, phi_deg);

    [lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);

    % actual transmitted pulse-1 uses pk1_exec
    [lfm_h1_exec, lfm_v1_exec, ~] = tx_polarize(lfm_temp, pk1_exec);
    [lfm_scatter1_exec, ~, ~] = build_target_echo(S, lfm_h1_exec, lfm_v1_exec, N_PRI, delay_samples);

    smsp1_exec = gen_smsp_polar(lfm_temp, pk1_exec, pj, T, B, fs, fc, M_smsp, ...
        delta_smsp, N_PRI, delay_samples, Pi);

    rsig1 = lfm_scatter1_exec + smsp1_exec;

    Ytr1 = rsig1(:, train_idx);
    phat = get_max_ev(Ytr1);
    phat = phat / max(norm(phat), 1e-12);

    %% ---------- Pulse 2: ideal design, mismatched execution ----------
    opts.opt2_pj = phat;
    pk2_ideal = gen_pk(polar_num, track2, S, opts);
    pk2_exec  = apply_phase_mismatch(pk2_ideal, phi_deg);

    % ideal pulse-2 target response used in receiver constraint
    [lfm_h2_id, lfm_v2_id, ~] = tx_polarize(lfm_temp, pk2_ideal);
    [lfm_scatter2_id, ~, ~] = build_target_echo(S, lfm_h2_id, lfm_v2_id, N_PRI, delay_samples);

    t = sum(lfm_scatter2_id(:, win), 2);
    if norm(t) == 0
        t = [1; 0];
    end

    % actual pulse-2 transmitted signal uses pk2_exec
    [lfm_h2_exec, lfm_v2_exec, ~] = tx_polarize(lfm_temp, pk2_exec);
    [lfm_scatter2_exec, ~, ~] = build_target_echo(S, lfm_h2_exec, lfm_v2_exec, N_PRI, delay_samples);

    smsp2_exec = gen_smsp_polar(lfm_temp, pk2_exec, pj, T, B, fs, fc, M_smsp, ...
        delta_smsp, N_PRI, delay_samples, Pi);

    %% ---------- Monte-Carlo over noise ----------
    s = zeros(1, MC);
    g = zeros(1, MC);

    for ii = 1:MC

        noise = [sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                 sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];

        rsig2 = lfm_scatter2_exec + smsp2_exec + noise;

        sinr_pre = get_sinr(rsig2, lfm_scatter2_exec);

        Ytr2  = rsig2(:, train_idx);
        Rhat  = (Ytr2 * Ytr2') / max(1, size(Ytr2, 2));
        Rhat  = (Rhat + Rhat') / 2;

        lambda = lambda_scale * trace(Rhat) / 2;
        Rdl    = Rhat + lambda * eye(2);

        C = [t, phat];
        f = [1; 0];

        Gm = C' * (Rdl \ C);
        if rcond(Gm) < 1e-6
            Gm = Gm + 1e-6 * eye(2);
        end

        w_lcmv = (Rdl \ C) * (Gm \ f);

        s(ii) = get_sinr(w_lcmv' * rsig2, w_lcmv' * lfm_scatter2_exec);
        g(ii) = s(ii) - sinr_pre;
    end

    sinr_avg = mean(s);
    gain_avg = mean(g);
end

%% =========================================================
function pk_out = apply_phase_mismatch(pk_in, phi_deg)
% Per-state relative H/V phase mismatch:
% actual Jones vector = [p_H; p_V * exp(j delta_i)], delta_i ~ U[-phi, phi]

    pk_out = pk_in;

    if phi_deg <= 0
        pk_out = norm_cols_local(pk_out);
        return;
    end

    phi_rad = phi_deg * pi / 180;
    M = size(pk_in, 2);

    delta = (2 * rand(1, M) - 1) * phi_rad;

    pk_out(2, :) = pk_out(2, :) .* exp(1j * delta);

    pk_out = norm_cols_local(pk_out);
end

%% =========================================================
function X = norm_cols_local(X)
    n = sqrt(sum(abs(X).^2, 1));
    n(n < 1e-12) = 1;
    X = bsxfun(@rdivide, X, n);
end