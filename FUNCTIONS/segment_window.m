% ===================== segment_window.m  =====================
function [seg_idx_train, seg_cnt] = segment_window(train_idx, seg_len_scalar, polar_num)
    
    delay_samples = train_idx(1) - 1;
    len = repmat(seg_len_scalar, 1, polar_num);
    seg_start = [1, cumsum(len(1:end - 1)) + 1];
    seg_end   = cumsum(len);
    
    K = polar_num;
    seg_idx_train = cell(1, K);
    seg_cnt       = zeros(1, K);
    for k = 1:K
        idxk = seg_start(k):seg_end(k);
        idxk = idxk + delay_samples;
        idxk = idxk(idxk >= train_idx(1) & idxk <= train_idx(end));
        seg_idx_train{k} = idxk;
        seg_cnt(k) = numel(idxk);
    end

end
