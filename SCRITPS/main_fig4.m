clc; clear; close all;

%%
config;
configplot;
set(groot, 'defaultScatterMarkerFaceAlpha', 0.1);

vec = S(2,:);
S = [1,                         vec(1) * exp(1j * vec(3));
     vec(1) * exp(1j * vec(3)), vec(2) * exp(1j * vec(4))];

%%
j_gamma = ((rand(1, 1) * 0.5) + 0.5) * (pi / 2);
j_eta   = (rand(1, 1) - 0.5) * (2 * pi);
j_gamma = deg2rad(48);
j_eta = deg2rad(120);
% j_gamma = deg2rad(82);
% j_eta = deg2rad(-60);
pj = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];

JSR_dB    = -5:5:5;
polar_num = 20;
train_len = 500;
MC        = 200;

tracks    = [2, 3, 13];
names     = {'RND','HV','RCL','GND','EQR', ...
    'TRG','ICO','FIB','TAM','COV', ...
    'OUR','RBT','OUR'};

train_idx = start_idx : (start_idx + train_len - 1);

[lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);

color_set = lines(numel(tracks));
res = struct();
for ij = 1:length(JSR_dB)
    jsr = JSR_dB(ij);
    Pi  = sigma * 10 ^ ((SNR + jsr) / 10);
    for it = 1:numel(tracks)
        tr = tracks(it);
        
        opts = struct('strategy','autotune','eta_floor',0.34,'deltaT',0.10, ...
            'rx_supp_frac',0.10,'jitter_deg',5,'phase_mode','uniform', ...
            'backproject', false);
        pk = build_pk(polar_num, tr, S, opts);

        [txH, txV]   = tx_polarize(lfm_temp, pk);
        [echo, ~, ~] = build_target_echo(S, txH, txV, N_PRI, delay_samples);

        ghat = zeros(MC, 1); 
        ehat = zeros(MC, 1);
        for m = 1:MC
            smsp = gen_smsp_polar(lfm_temp, pk, pj, T, B, fs, fc, ...
                M_smsp, 0, N_PRI, delay_samples, Pi);
            
            noise = [sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                     sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];
        
            rsig = echo + smsp + noise;
                 
            Ytr = rsig(:, train_idx);
            phat = get_max_ev(Ytr);

            [g, e] = vec2angles(phat);
            ghat(m) = g; 
            ehat(m) = e;
        end
        
        res(it).name  = names{tr};
        res(it).gamma = ghat; 
        res(it).eta   = ehat;
    end

    figure('Color','w'); hold on; box on; grid on;
    
    hTrue = plot(rad2deg(j_gamma), rad2deg(wrapToPi(j_eta)), ...
        'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
        'DisplayName', 'TRUE');
    
    hLeg = gobjects(numel(tracks), 1);
    
    for it = 1:numel(tracks)
        
        scatter(rad2deg(res(it).gamma), rad2deg(wrapToPi(res(it).eta)), ...
            14, 'filled', 'MarkerFaceAlpha', 0.25, 'MarkerEdgeAlpha', 0.2, ...
            'MarkerFaceColor', color_set(it, :), 'HandleVisibility', 'off');
        
        hLeg(it) = plot(nan,nan, 'o', ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', color_set(it,:), ...
            'MarkerEdgeColor', color_set(it,:), ...
            'LineStyle', 'none', ...
            'DisplayName', names{tracks(it)});
    end
    
    xlabel('\gamma (deg)', 'FontName', 'Times New Roman');
    ylabel('\eta (deg)',   'FontName', 'Times New Roman');
    xlim([0 90]); ylim([-180 180]); yticks(-180:60:180);
    
    legend([hTrue; hLeg], [{'TRUE'}, names(tracks)], 'Location','southwest');

end


%%
function [g, e] = vec2angles(p)
    p = p(:); 
    p = p / max(norm(p), 1e-12); 
    p = p * exp(-1j * angle(p(1)));
    g = atan2(abs(p(2)), abs(p(1)));
    e = angle(p(2));
end
