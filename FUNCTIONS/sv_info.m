function [S, Sigma, max_cs_bound] = sv_info(k, r, phx, phv)
% SV_INFO  构造 2x2 极化散射矩阵 S，并给出奇异值与余弦相似度上界
%   S = [1,  k*exp(1j*phx);
%        k*exp(1j*phx),  r*exp(1j*phv)];
%
% 输出:
%   S            : 2x2 复矩阵
%   svals        : 奇异值 [sigma1; sigma2]（降序）
%   max_cs_bound : 最大余弦相似度的理论上界  sigma1 / ||S||_F

    S = [1,                  k * exp(1j * phx);
         k * exp(1j * phx),  r * exp(1j * phv)];

    a = 1 + k ^ 2;
    b = k * exp(1j * phx) + k * r * exp(1j * (phv - phx));
    c = k * exp(1j * -phx) + k * r * exp(1j * (phx - phv));
    d = k ^ 2 + r ^ 2;
    
    t1 = (a + d) / 2;
    t2 = sqrt((a + d) ^ 2 / 4 - (a * d - b * c));
    lambda_1 = t1 + t2;
    lambda_2 = t1 - t2;
    
    sigma1 = sqrt(lambda_1);
    sigma2 = sqrt(lambda_2);
    
    Sigma = [sigma1, sigma2];
    
    max_cs_bound = sigma1 / norm(S, 'fro');
    
end
