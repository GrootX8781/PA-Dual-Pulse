clc; clear; close all;

%%
config;
configplot;

load psm.mat;

vec = s1(2,:);
k   = vec(1);  rs = vec(2);  phx = vec(3);  phv = vec(4);
S   = [1,                 k  * exp(1j * phx);
       k * exp(1j * phx), rs * exp(1j * phv)];

%%
j_gamma = pi / 3;
j_eta   = pi / 2;
% j_gamma = ((rand(1, 1) * 0.5) + 0.5) * (pi / 2);
% j_eta   = (rand(1, 1) - 0.5) * (2 * pi);
pj      = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];

polar_num = 30;
train_len = N_temp;

JSR = -5:5:30;
MC  = 100;

%%
track1s = [14, 2, 3, 4, 10];
track2s = [14, 2, 3, 4, 11];

% track1s = [12, 13, 2, 3, 4, 10];
% track2s = [12, 13, 2, 3, 4, 11];

sinr = zeros(length(track1s), length(JSR));
gain = zeros(length(track1s), length(JSR));

for it = 1:length(track1s)
    
    fprintf('TESTING %2d /%2d\n', it, length(track1s));
    
    track1 = track1s(it);
    track2 = track2s(it);
    
    [s, g] = run_once2(JSR, polar_num, track1, track2, S, ...
        train_len, j_gamma, j_eta, MC);
    
    sinr(it, :) = s;
    gain(it, :) = g;
    
end

%%
names = {'RND', 'HV', 'RCL', 'GND', 'EQR', ...
    'TRG', 'ICO', 'FIB', 'COV', 'OPT1', ...
    'OPT2', 'H极化', 'V极化', 'LIN'};

figure('Color', 'w'); hold on; grid on;
for k = 1:numel(track1s)
    name1 = names{track1s(k)};
    name2 = names{track2s(k)};
    if k~= numel(track1s)
        name = name1;
    else
        name = 'OUR';
    end
    plot(JSR, sinr(k, :), 'LineWidth', lw, ...
        'Color', cols_ln{k}, 'LineStyle', ls{k}, ...
        'DisplayName', name);
end
xlabel('JSR (dB)', 'FontName', 'Times New Roman'); 
ylabel('SINR (dB)', 'FontName', 'Times New Roman');
legend('Location', 'southeast');
ylim([-25 5])
% ylim([round(min(sinr(:)) - 3) round(max(sinr(:)) + 3)])
set(gcf, 'Position', [300 300 800 600]);

figure('Color', 'w'); hold on; grid on;
for k = 1:numel(track1s)
    name1 = names{track1s(k)};
    name2 = names{track2s(k)};
    if k ~= numel(track1s)
        name = name1;
    else
        name = 'OUR';
    end
    plot(JSR, gain(k, :), 'LineWidth', lw, ...
        'Color', cols_ln{k}, 'LineStyle', ls{k}, ...
        'DisplayName', name);
end
xlabel('JSR (dB)', 'FontName', 'Times New Roman'); 
ylabel('SINR Gain(dB)', 'FontName', 'Times New Roman');
legend('Location', 'southeast');
ylim([-15 30])
set(gcf, 'Position', [200 200 800 600]);

%%
function [sinr, gain] = run_once2(JSR_list, polar_num, track1, track2, S, ...
    train_len, j_gamma, j_eta, MC)

    config;
    
    SNR      = 10;
    Ps       = 1 * 10 .^ (SNR / 10);
    
    pj = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];
    
    sinr = zeros(1, length(JSR_list));
    gain = zeros(1, length(JSR_list));
    lambda_scale = 5e-2;
       
    win       = start_idx : end_idx;
    train_idx = start_idx : (start_idx + train_len - 1);
    
    for ij = 1:length(JSR_list)
        
        jsr = JSR_list(ij);
        Pi = 10 .^ ((jsr + SNR) / 10);
        
        pk1 = gen_pk(polar_num, track1, S);
        [lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);
        [lfm_h, lfm_v, ~] = tx_polarize(lfm_temp, pk1);
        [lfm_scatter1, ~, ~] = build_target_echo(S, lfm_h, lfm_v, N_PRI, delay_samples);
        
        smsp1 = gen_smsp_polar(lfm_temp, pk1, pj, T, B, fs, fc, M_smsp, ...
            delta_smsp, N_PRI, delay_samples, Pi);
        
        rsig1 = lfm_scatter1 + smsp1;

        Ytr = rsig1(:, train_idx);
        phat = get_max_ev(Ytr);
        
        opts.opt2_pj = phat;
        pk2 = gen_pk(polar_num, track2, S, opts);

        [lfm_h, lfm_v, ~] = tx_polarize(lfm_temp, pk2);
        [lfm_scatter2, ~, ~] = build_target_echo(S, lfm_h, lfm_v, N_PRI, delay_samples);
        t = sum(lfm_scatter2(:, win), 2); if norm(t)==0, t=[1;0]; end
        
        smsp2 = gen_smsp_polar(lfm_temp, pk2, pj, T, B, fs, fc, M_smsp, ...
            delta_smsp, N_PRI, delay_samples, Pi);

        s = zeros(1, MC);
        g = zeros(1, MC);
        for ii = 1:MC
            noise = [sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                     sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];
                 
            rsig2 = lfm_scatter2 + smsp2 + noise;
            
            sinr_pre = get_sinr(rsig2, lfm_scatter2);
            
            Ytr    = rsig2(:, train_idx);
            Rhat   = (Ytr * Ytr') / max(1, size(Ytr, 2)); 
            Rhat   = (Rhat + Rhat') / 2;
            lambda = lambda_scale * trace(Rhat) / 2;
            Rdl    = Rhat + lambda * eye(2);

            C = [t, phat]; 
            f = [1; 0];
            Gm = C' * (Rdl \ C); if rcond(Gm)<1e-6, Gm = Gm + 1e-6*eye(2); end
            w_lcmv = (Rdl \ C) * (Gm \ f);
            s(ii) = get_sinr(w_lcmv' * rsig2, w_lcmv' * lfm_scatter2);
            g(ii) = s(ii) - sinr_pre;
        end
        
        sinr(ij) = mean(s);
        gain(ij) = mean(g);

    end
end

