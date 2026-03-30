% ===================== gen_lfm.m  =====================
function [lfm_temp, N_temp] = gen_lfm(T, B, fs, fc, Ps)
% 生成一次脉冲的 LFM 基带（含载频），幅度按 Ps 归一
    K = B/T;
    N_temp = round(T*fs);
    t = (0:N_temp-1)/fs;
%     t = linspace(-T / 2, T / 2, N_temp);
    lfm_temp = sqrt(Ps) * exp(1j*pi*K*t.^2) .* exp(1j*2*pi*fc*t);
end
