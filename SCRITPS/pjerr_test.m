clc; clear; close all

%%
config;
configplot;

load psm.mat;

s = s1(2, :);
[s1, s2, s3, s4] = deal(s(1), s(2), s(3), s(4));
[~, ~, bound] = sv_info(s1, s2, s3, s4);
S = [1,                 s1 * exp(1j * s3);
     s1 * exp(1j * s3), s2 * exp(1j * s4)];

Ng = 91;
Ne = 91;
gamma_bin = linspace(pi / 4, pi / 2, Ng);
eta_bin   = linspace(-pi, pi, Ne);

MC = 1;

%%
polar_num = 20;
K         = 2000;
JSR       = -5:5:30;
names = {'RND', 'HV', 'RCL', 'GND', 'EQR', ...
    'TRG', 'ICO', 'FIB', 'COV', 'OUR', ...
    'OPT2', 'single-H', 'single-V', 'LIN'};

test_type = 1;
test_case = [14, 2, 3, 4, 10];
% test_case = [5 11];
err = zeros(length(test_case), length(JSR));

fprintf('CASE PARAMETER TEST\n');
fprintf('NUM    NAME   lamminT   kappaT     M1      M2     lamminR   kappaR\n');
for it = 1:length(test_case)
    
    track = test_case(it);
    pk = gen_pk(polar_num, track, S);
    info = pk_info(pk, S);
    fprintf('%2d %8s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
        it, names{track}, info.lamminT, info.kappaT, info.M1, ...
        info.M2, info.lamminR, info.kappaR);
    
    err(it, :) = run_once(1, JSR, polar_num, track, ...
        S, K, gamma_bin, eta_bin, MC);
    
end


%%
function err = run_once(test_type, jsr_list, polar_num, track, ...
    S, train_len, gamma_bin, eta_bin, MC)

   config;

   pk = gen_pk(polar_num, track, S);
   [lfm_temp, ~] = gen_lfm(T, B, fs, fc, Ps);
   [lfm_h, lfm_v, ~] = tx_polarize(lfm_temp, pk);
   [lfm_scatter, ~, ~] = build_target_echo(S, lfm_h, lfm_v, N_PRI, delay_samples);
   
    seg_len_scalar     = round(N_temp / polar_num);
    train_idx          = start_idx : (start_idx + train_len - 1);
    [seg_idx_train, ~] = segment_window(train_idx, seg_len_scalar, polar_num);

    err = zeros(size(jsr_list));

    for ij = 1:length(jsr_list)
        jsr     = jsr_list(ij);
        Pi_smsp = 10 ^ ((jsr + SNR) / 10);

        acc = zeros(1, MC * numel(gamma_bin) * numel(eta_bin));
        ptr = 0;

        for mc = 1:MC
            
            noise = [sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI));
                     sqrt(1 / 2) * (randn(1, N_PRI) + 1j * randn(1, N_PRI))];

            for ig = 1:numel(gamma_bin)
                j_gamma = gamma_bin(ig);
                for ie = 1:numel(eta_bin)
                    ptr   = ptr + 1;
                    j_eta = eta_bin(ie);

                    pj = [cos(j_gamma); sin(j_gamma) * exp(1j * j_eta)];
                    pj = pj / max(norm(pj),1e-12); 
                    pj = pj * exp(-1j * angle(pj(1)));

                    smsp = gen_smsp_polar(lfm_temp, pk, pj, T, B, fs, fc, M_smsp, ...
                        delta_smsp, N_PRI, delay_samples, Pi_smsp);

%                     rsig = lfm_scatter + smsp + noise;
                    rsig = lfm_scatter + smsp;
                    
                    if test_type == 1
                        Ytr = rsig(:, train_idx);
                        phat = get_max_ev(Ytr);
                    elseif test_type == 2
                        opts = struct('lfm_temp',lfm_temp,'start_idx',start_idx, ...
                                      'iters_refine',1,'q_bell',1.0,'debug',false);
                        [phat, ~] = pj_est_test(rsig, pk, [], seg_idx_train, train_idx, [], [], [], opts);
                    end
                    er = 1 - (abs(phat' * pj) ^ 2) / (norm(phat)^2 * norm(pj)^2);
                    acc(ptr) = max(0, real(er));
                end
            end
        end
        err(ij) = exp(mean(log(acc + 1e-15)));
    end
end



