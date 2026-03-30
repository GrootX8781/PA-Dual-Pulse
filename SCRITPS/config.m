% ===================== config.m  =====================
% 统一参数配置（按需修改）
c   = 3.0e8;
T   = 10e-6;
B   = 50e6;
fs  = 4 * B;
fc  = 15e9;
PRI = 100e-6;
lambda = c / fc;
d = lambda / 2;
d_lambda = d / lambda;

r        = 9000;           % 目标距离
SNR      = 0;             % 目标 SNR（定义在发射模板上）
JSR      = 30;             % 干扰相对 SNR（dB）
sigma    = 1;
Ps       = sigma * 10 ^ (SNR / 10);
Pi_smsp  = sigma * 10 ^ ((SNR + JSR) / 10);

% DRFM/SMSP
M_smsp     = 3;         % 压缩倍数
delta_smsp = 0;

Ts = 2e-6;
tau = 0.5e-6;
cycle = tau / Ts;
INR_isrj = SNR + JSR;
delta_isrj = 0;


delay_samples = round(2 * r / c * fs);
N_PRI         = round(PRI * fs);
N_temp        = round(T * fs);
start_idx     = delay_samples + 1;
end_idx       = start_idx + N_temp - 1;


t     = (0:(N_PRI - 1)) / fs;
range = t * c / 2;

