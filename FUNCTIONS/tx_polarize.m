% ===================== tx_polarize.m  =====================
function [lfm_h, lfm_v, seg_len] = tx_polarize(lfm_temp, pk)
% 把 LFM 模板按分段极化映射到 H/V 两路
    N_temp  = numel(lfm_temp);
    K       = size(pk,2);
    seg_len = round(N_temp / K);

    lfm_h = zeros(1, N_temp); 
    lfm_v = zeros(1, N_temp);
    for i = 1:K
        a = (i-1)*seg_len + 1;
        if i < K, b = a + seg_len - 1; else, b = N_temp; end
        lfm_h(a:b) = pk(1,i) * lfm_temp(a:b);
        lfm_v(a:b) = pk(2,i) * lfm_temp(a:b);
    end
end
