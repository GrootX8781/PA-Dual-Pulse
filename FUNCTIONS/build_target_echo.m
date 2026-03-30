% ===================== build_target_echo.m  =====================
function [lfm_scatter, start_idx, end_idx] = build_target_echo(S, lfm_h, lfm_v, N_PRI, delay_samples)
% 把两路发射模板按时延嵌入 PRI 轴，再乘极化散射矩阵 S
    N_temp = numel(lfm_h);
    start_idx = delay_samples + 1;
    end_idx   = start_idx + N_temp - 1;

    lfm_echo = zeros(2, N_PRI);
    lfm_echo(1, start_idx:end_idx) = lfm_h;
    lfm_echo(2, start_idx:end_idx) = lfm_v;

    lfm_scatter = S * lfm_echo;   % 目标极化散射
end
