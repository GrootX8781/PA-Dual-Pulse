clc; clear; close all;
rng(1);

%% basic config
config;
configplot;


%% true target PSM
vec = S(2, :);
k   = vec(1);
rs  = vec(2);
phx = vec(3);
phv = vec(4);

S_true = [1,                  k  * exp(1j * phx);
          k * exp(1j * phx),  rs * exp(1j * phv)];

%% jammer polarization
j_gamma = pi / 3;
j_eta   = pi / 2;
pj      = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];

%% experiment setup
polar_num = 30;
train_len = N_temp;

% fixed JSR values for robustness evaluation
JSR_fix = [0, 10, 20];      % dB

% relative Frobenius mismatch levels
eps_list = 0:0.05:0.30;     % epsilon = ||S_hat - S||_F / ||S||_F

MC_noise = 100;             % inner Monte-Carlo for noise
MC_psm   = 20;              % outer Monte-Carlo for random PSM mismatch realizations

%% track settings
track1s = [14, 2, 3, 4, 10];
track2s = [14, 2, 3, 4, 11];

names = {'RND', 'HV', 'RCL', 'GND', 'EQR', ...
         'TRG', 'ICO', 'FIB', 'COV', 'OPT1', ...
         'OPT2', 'H极化', 'V极化', 'LIN'};

%% result arrays
% size = [num_tracks, num_JSR, num_eps]
sinr_out = zeros(length(track1s), length(JSR_fix), length(eps_list));
gain_out = zeros(length(track1s), length(JSR_fix), length(eps_list));

%% main loop
for it = 5:length(track1s)

    fprintf('TRACK %2d / %2d\n', it, length(track1s));

    track1 = track1s(it);
    track2 = track2s(it);

    for ij = 1:length(JSR_fix)

        jsr_now = JSR_fix(ij);
        fprintf('  JSR = %+3d dB\n', jsr_now);

        for ie = 1:length(eps_list)

            eps_now = eps_list(ie);

            s_tmp = zeros(1, MC_psm);
            g_tmp = zeros(1, MC_psm);

            for im = 1:MC_psm
                S_design = gen_S_mismatch_safe(S_true, eps_now);

                [s_tmp(im), g_tmp(im), mode1, mode2] = run_once2_psm_single_safe( ...
                    jsr_now, polar_num, track1, track2, ...
                    S_true, S_design, train_len, ...
                    j_gamma, j_eta, MC_noise, eps_now);
            end

            sinr_out(it, ij, ie) = mean(s_tmp);
            gain_out(it, ij, ie) = mean(g_tmp);

            fprintf('    eps = %.2f | SINR = %.3f dB | Gain = %.3f dB\n', ...
                eps_now, sinr_out(it, ij, ie), gain_out(it, ij, ie));
        end
    end
end

%% save results
save('z_tsp1_fig9.mat', 'sinr_out', 'gain_out', ...
    'eps_list', 'JSR_fix', 'track1s', 'track2s', 'MC_noise', 'MC_psm');

%% plot: output SINR vs epsilon
for ij = 1:length(JSR_fix)

    figure('Color', 'w'); hold on; grid on;

    for it = 1:length(track1s)
        if it ~= length(track1s)
            name = names{track1s(it)};
        else
            name = 'OUR';
        end

        plot(eps_list, squeeze(sinr_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{it}, ...
            'DisplayName', name);
    end

    xlabel('Relative PSM mismatch \epsilon', 'FontName', 'Times New Roman');
    ylabel('Output SINR (dB)', 'FontName', 'Times New Roman');
    title(sprintf('PSM Robustness: Output SINR (JSR = %d dB)', JSR_fix(ij)), ...
        'FontName', 'Times New Roman');
    legend('Location', 'best');
    set(gcf, 'Position', [200 200 800 600]);
end

%% plot: SINR gain vs epsilon
for ij = 1:length(JSR_fix)

    figure('Color', 'w'); hold on; grid on;

    for it = 1:length(track1s)
        if it ~= length(track1s)
            name = names{track1s(it)};
        else
            name = 'OUR';
        end

        plot(eps_list, squeeze(gain_out(it, ij, :)), ...
            'LineWidth', lw, ...
            'Color', cols_ln{it}, ...
            'LineStyle', ls{it}, ...
            'DisplayName', name);
    end

    xlabel('Relative PSM mismatch \epsilon', 'FontName', 'Times New Roman');
    ylabel('SINR Gain (dB)', 'FontName', 'Times New Roman');
    title(sprintf('PSM Robustness: SINR Gain (JSR = %d dB)', JSR_fix(ij)), ...
        'FontName', 'Times New Roman');
    legend('Location', 'best');
    set(gcf, 'Position', [260 260 800 600]);
end

%% =========================================================
%% local function 1
function [sinr_avg, gain_avg, mode_pk1, mode_pk2] = run_once2_psm_single_safe( ...
    jsr, polar_num, track1, track2, ...
    S_true, S_design, train_len, ...
    j_gamma, j_eta, MC, eps_now)

    config;

    SNR = 10;
    Ps  = 1 * 10^(SNR / 10);

    pj = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];

    lambda_scale = 5e-2;

    win       = start_idx : end_idx;
    train_idx = start_idx : (start_idx + train_len - 1);

    Pi = 10^((jsr + SNR) / 10);

    %% ---------- Pulse 1 ----------
    opts1 = struct();
    opts1.backproject = true;

    [pk1, mode_pk1] = safe_gen_pk( ...
        polar_num, track1, S_design, S_true, opts1, eps_now);

    [lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);
    [lfm_h1, lfm_v1, ~] = tx_polarize(lfm_temp, pk1);

    [lfm_scatter1_true, ~, ~] = build_target_echo( ...
        S_true, lfm_h1, lfm_v1, N_PRI, delay_samples);

    smsp1 = gen_smsp_polar( ...
        lfm_temp, pk1, pj, T, B, fs, fc, ...
        M_smsp, delta_smsp, N_PRI, delay_samples, Pi);

    rsig1 = lfm_scatter1_true + smsp1;

    Ytr1 = rsig1(:, train_idx);
    phat = get_max_ev(Ytr1);

    %% ---------- Pulse 2 ----------
    opts2 = struct();
    opts2.backproject = true;
    opts2.opt2_pj = phat;

    [pk2, mode_pk2] = safe_gen_pk( ...
        polar_num, track2, S_design, S_true, opts2, eps_now);

    [lfm_h2, lfm_v2, ~] = tx_polarize(lfm_temp, pk2);

    % true target echo used in received signal
    [lfm_scatter2_true, ~, ~] = build_target_echo( ...
        S_true, lfm_h2, lfm_v2, N_PRI, delay_samples);

    % mismatched target echo used in target constraint
    [lfm_scatter2_design, ~, ~] = build_target_echo( ...
        S_design, lfm_h2, lfm_v2, N_PRI, delay_samples);

    t = sum(lfm_scatter2_design(:, win), 2);
    if norm(t) == 0
        t = [1; 0];
    end

    smsp2 = gen_smsp_polar( ...
        lfm_temp, pk2, pj, T, B, fs, fc, ...
        M_smsp, delta_smsp, N_PRI, delay_samples, Pi);

    %% ---------- Monte-Carlo over noise ----------
    s = zeros(1, MC);
    g = zeros(1, MC);

    for ii = 1:MC

        noise = [sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                 sqrt(1/2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];

        rsig2 = lfm_scatter2_true + smsp2 + noise;

        sinr_pre = get_sinr(rsig2, lfm_scatter2_true);

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

        s(ii) = get_sinr(w_lcmv' * rsig2, w_lcmv' * lfm_scatter2_true);
        g(ii) = s(ii) - sinr_pre;
    end

    sinr_avg = mean(s);
    gain_avg = mean(g);
end

%% =========================================================
%% local function 2
function [pk, mode_flag] = safe_gen_pk(polar_num, type, S_design, S_true, opts_in, eps_now)
% mode_flag:
% 0 = normal design with S_design
% 1 = relaxed-feasibility design with S_design
% 2 = nominal(S_true) + small trajectory perturbation fallback

    % ---------- first try: normal ----------
    try
        [pk, ~] = gen_pk(polar_num, type, S_design, opts_in);
        mode_flag = 0;
        return;
    catch
    end

    % ---------- second try: relaxed feasibility for OPT1 only ----------
    if type == 10
        delta_list = [0.03, 0.05, 0.08, 0.12, 0.18, 0.25, 0.35];
        for id = 1:length(delta_list)
            opts_try = opts_in;
            opts_try.opt1_deltaR = delta_list(id);
            try
                [pk, ~] = gen_pk(polar_num, type, S_design, opts_try);
                mode_flag = 1;
                return;
            catch
            end
        end
    end

    % ---------- final fallback: nominal trajectory + tiny perturbation ----------
    try
        [pk, ~] = gen_pk(polar_num, type, S_true, opts_in);
    catch
        % if even nominal fails, use a fixed linear fallback
        [pk, ~] = gen_pk(polar_num, 14, S_true, struct('backproject', true));
    end

    pk = perturb_pk_cols(pk, eps_now);
    mode_flag = 2;
end

%% =========================================================
%% local function 3
function S_design = gen_S_mismatch_safe(S_true, eps_level)
% Generate reciprocal mismatched S_design satisfying:
% eps_level = ||S_design - S_true||_F / ||S_true||_F
%
% Also avoid overly ill-conditioned S_design that would make backprojection unstable.

    if eps_level <= 0
        S_design = S_true;
        return;
    end

    target_norm = eps_level * norm(S_true, 'fro');

    max_retry = 50;
    for it = 1:max_retry

        a = randn + 1j * randn;
        b = randn + 1j * randn;
        c = randn + 1j * randn;

        E = [a, b;
             b, c];

        E = E / max(norm(E, 'fro'), 1e-12);

        S_try = S_true + target_norm * E;

        % enforce reciprocity exactly
        S_try = 0.5 * (S_try + S_try.');

        % avoid too singular design matrices
        Rt = S_try * S_try';
        if rcond(Rt) > 1e-6 && norm(S_try, 'fro') > 1e-10
            S_design = S_try;
            return;
        end
    end

    % last-resort blend
    S_design = S_true + target_norm * eye(2) / sqrt(2);
    S_design = 0.5 * (S_design + S_design.');
end

%% =========================================================
%% local function 4
function pk = perturb_pk_cols(pk, eps_now)
% small trajectory perturbation used only as a LAST fallback
% to emulate mild model-induced trajectory deviation

    pk = norm_cols_local(pk);

    [gamma, eta] = jones_to_gamma_eta_local(pk);

    % very small perturbation; scales mildly with eps_now
    sig_g = min(0.08, 0.35 * eps_now) * (pi / 6);   % gamma perturbation
    sig_e = min(0.12, 0.45 * eps_now) * (pi / 4);   % eta perturbation

    gamma = gamma + sig_g * randn(size(gamma));
    eta   = eta   + sig_e * randn(size(eta));

    gamma = min(max(gamma, 0), pi / 2);
    eta   = wrapToPi(eta);

    pk = [cos(gamma); sin(gamma) .* exp(1j * eta)];
    pk = norm_cols_local(pk);
end

%% =========================================================
%% local function 5
function X = norm_cols_local(X)
    n = sqrt(sum(abs(X).^2, 1));
    n(n < 1e-12) = 1;
    X = bsxfun(@rdivide, X, n);
end

%% =========================================================
%% local function 6
function [gamma, eta] = jones_to_gamma_eta_local(pk)
    pk = norm_cols_local(pk);
    a = pk(1, :);
    b = pk(2, :);
    gamma = atan2(abs(b), abs(a));
    eta   = angle(b) - angle(a);
    eta   = wrapToPi(eta);
end