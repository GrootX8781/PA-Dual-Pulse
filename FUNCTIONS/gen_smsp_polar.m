% ===================== gen_smsp_polar.m  =====================
function smsp = gen_smsp_polar(lfm_temp, pk, pj, T, B, fs, fc, ...
                               M_smsp, delta_smsp, N_PRI, delay_samples, Pi_smsp)
% 极化 SMSP：pj.' 接收投影 → 复制压缩 → pj 发射；最终功率定标到 Pi_smsp
    c = 3e8; %#ok<NASGU>  % 仅用于 idx 计算时沿用主脚本常量时可移除
    T_sub = T / M_smsp;
    K_sub = B / T_sub;
    N_sub = round(T_sub*fs);
    t_sub = (0:N_sub-1)/fs;

    % 单位功率的压缩副本
    sub_temp = exp(1j*pi*K_sub*t_sub.^2 + 1j*2*pi*fc*t_sub);

    % 构造 2×N_sub 的分段极化模板（沿用 pk）
    K = size(pk,2);
    nn_sub   = round(N_sub / K);
    sub_base = zeros(2, N_sub);
    for i = 1:K
        a = (i-1)*nn_sub + 1; b = min(a+nn_sub-1, N_sub);
        sub_base(1, a:b) = pk(1,i) * sub_temp(a:b);
        sub_base(2, a:b) = pk(2,i) * sub_temp(a:b);
    end

    % 接收：pj.' * sub_base  （模型要求：转置、非共轭）
    polar_sub = pj.' * sub_base;               % 1×N_sub

    % 复制压缩并裁剪到 N_temp
    temp = repmat(polar_sub, 1, M_smsp);       % 1×(M*N_sub)
    N_temp = numel(lfm_temp);
    temp = temp(:, 1:N_temp);

    % 发射：pj * temp → 2×N_temp
    temp = pj * temp;

    % 归一化到单位平均功率，再按 Pi_smsp 定标
    P0 = mean(sum(abs(temp).^2, 1)); 
    if P0 > 0, temp = temp / sqrt(P0); end
    temp = sqrt(Pi_smsp) * temp;

    % 放入 PRI 轴（相对目标时延 + 额外偏移）
    idx0    = delay_samples + round(2*delta_smsp/3e8 * fs);
    j_start = idx0 + 1; 
    j_end   = j_start + N_temp - 1;

    smsp = zeros(2, N_PRI);
    if j_start>=1 && j_end<=N_PRI
        smsp(:, j_start:j_end) = temp;
    else
        % 边界安全写入
        a = max(1, j_start); b = min(N_PRI, j_end);
        if a<=b
            sa = 1 + (a - j_start); sb = sa + (b - a);
            smsp(:, a:b) = temp(:, sa:sb);
        end
    end
end
