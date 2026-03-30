function sinr = get_sinr(sig, tar)
    e_s  = sum(abs(tar) .^ 2, 'all');
    e_jn = sum(abs(sig - tar) .^ 2, 'all');
    sinr = 10 * log10(e_s / max(e_jn, eps));
end