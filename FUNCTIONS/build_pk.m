function [pk, pk_rad] = build_pk(polar_num, polar_type, S, opts)
% 构造发射极化轨道 pk (2xK Jones 向量) 及其 (gamma, eta)
% polar_type:
%   1: 随机
%   2: H→V 线极化（eta=0，gamma 等分）
%   3: RC→(穿过45°线极化)→LC（gamma=pi/4，eta 从 +pi/2 到 -pi/2）
%   4: gamma 与 eta 同时等分扫描
%   5: Fibonacci 球均匀覆盖 → (可选) S? 反投影
%   6: 正二十面体顶点集 → (可选) S? 反投影   <-- 修复了等间隔取整索引
%   7: 等角环（equator ring, gamma=pi/4，eta 等分）→ (可选) S? 反投影
%   8: 双环 two-ring（gamma=pi/4±Δ，eta 等分）→ (可选) S? 反投影
%   9: Tammes解

    if nargin < 3 || isempty(S), S = []; end
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'rng_seed'),      opts.rng_seed = 1; end
    if ~isfield(opts,'backproject'),   opts.backproject = true; end
    if ~isfield(opts,'tworing_delta'), opts.tworing_delta = 15; end  % degree
    if ~isfield(opts,'conv'),          opts.conv = 'RHC+Z'; end

    K = polar_num;
    do_backproj = opts.backproject;

    switch polar_type
        case 1  % 随机
            rng(4);
            gamma  = sort(rand(1, K) * (pi/2));
            eta    = (rand(1, K) - 0.5) * 2 * pi;
            pk     = [cos(gamma); sin(gamma).*exp(1j*eta)];
            pk     = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 2  % H→V 线极化
            gamma  = linspace(0, pi/2, K);
            eta    = zeros(1, K);
            pk     = [cos(gamma); sin(gamma)];
            pk     = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 3  % RC→45°→LC
            gamma = (pi/4) * ones(1,K);
            eta   = linspace(+pi/2, -pi/2, K);
            pk    = [cos(gamma); sin(gamma).*exp(1j*eta)];
            pk    = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 4  % gamma 与 eta 同时等分扫描
            gamma = linspace(0, pi/2, K);
            eta   = linspace(-pi, pi, K);
            pk    = [cos(gamma); sin(gamma).*exp(1j*eta)];
            pk    = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
        case 5  % EQR 赤道环
            eta   = linspace(0, 2*pi, K+1); eta(end) = [];
            gamma = (pi/4) * ones(1, K);
            U = [cos(2*gamma); sin(2*gamma).*cos(eta); sin(2*gamma).*sin(eta)];
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 6  % TRG 双环
            d = opts.tworing_delta * pi/180;     % 自己转弧度，避免依赖 deg2rad
            d = 15 * pi / 180;
            d = 15*pi/180;
            phi = 0.7*(pi/ (K - floor(K/2)));    % 非对称相移
            K1 = floor(0.55*K); K2 = K-K1;       % 不等分
            eta1 = linspace(0, 2*pi, K1+1); eta1(end)=[];
            eta2 = linspace(phi, 2*pi+phi, K2+1); eta2(end)=[];
%             K1 = floor(K/2); K2 = K - K1;
%             eta1 = linspace(0, 2*pi, K1+1); eta1(end) = [];
%             eta2 = linspace(pi/K2, 2*pi+pi/K2, K2+1); eta2(end) = [];
            g1 = (pi/4 - d)*ones(1,K1); g2 = (pi/4 + d)*ones(1,K2);
            U1 = [cos(2*g1); sin(2*g1).*cos(eta1); sin(2*g1).*sin(eta1)];
            U2 = [cos(2*g2); sin(2*g2).*cos(eta2); sin(2*g2).*sin(eta2)];
            U  = [U1, U2];
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);

        case 7  % 正二十面体顶点 + (可选) S? 反投影  —— 修复索引
            U0 = icosahedron_vertices();         % 3xN (N≈20)
            n  = size(U0,2);
            if K <= n
                % 用整型等间隔索引，并去重，不足则补齐
                idx = round(linspace(1, n, K));
                idx(idx < 1) = 1; idx(idx > n) = n;
                idx = unique(idx, 'stable');
                if numel(idx) < K
                    pool = setdiff(1:n, idx, 'stable');
                    need = K - numel(idx);
                    idx  = [idx, pool(1:need)];
                end
                U = U0(:, idx);
            else
                rep = ceil(K/n);
                U = repmat(U0, 1, rep);
                U = U(:, 1:K);
            end
            U  = norm_cols(U);
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);

        case 8  % Fibonacci 球覆盖 + (可选) S? 反投影
            U = fibonacci_sphere(K);             % 3xK, (x,y,z)
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 9
            [U, ~, ~, ~] = tammes_points(K, ...
                'restarts', 25, 'display','off', 'seed', 42, 'use_table', false);
            pk = stokes3_to_jones_and_backproj(U.', S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 10
            load('best_cover.mat');
            if K < 4 || K > 30
                error('WRONG POLAR NUM');
            end
            U = cell2mat(best_cover(K - 3));
            pk = stokes3_to_jones_and_backproj(U.', S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 11  % OPT–Ring：κR≤δ（硬约束）、Tx 均匀覆盖、逐对相位使 Hbar=0
            if isempty(S), error('case 14 需要 S'); end
            if ~isfield(opts,'deltaR'), opts.deltaR = 0.02; end
            if ~isfield(opts,'oriN'),   opts.oriN   = max(800, 160*polar_num); end
            K = polar_num;
            if mod(K,2)~=0
                warning('K 需为偶数；已减 1。');  K = K-1;
            end

            % === 1) 先求 G*（κR≤δ 的最小 κT 解）：复用 case13 思路 ===
            [u_star, lam_star] = solve_u_lambda_minKappaT(S, opts.deltaR, opts.oriN);
            beta = 1 - 2*lam_star;                  % β ∈ [0,1]
            % u_star -> n（Stokes 方向）
            gu = atan2(abs(u_star(2)), abs(u_star(1)));
            et = angle(u_star(2)) - angle(u_star(1));
            n  = [cos(2*gu); sin(2*gu)*cos(et); sin(2*gu)*sin(et)];
            n  = n / max(norm(n),1e-12);

            % === 2) 轴向正交基（把 n 当作“z 轴”） ===
            % 选一条与 n 不共线的参考向量
            if abs(n(1) < 0.9)
                ref = [1; 0; 0];
            else
                ref = [0; 1; 0];
            end
            e1  = cross(n, ref); e1 = e1 / max(norm(e1),1e-12);
            e2  = cross(n, e1);  e2 = e2 / max(norm(e2),1e-12);

            % === 3) 环 + 中心 的解析构造 ===
            Q  = floor(K/4);         % 环上“方位对”数量（每对 4 列：两方位×两相位）
            R  = K - 4*Q;            % 0 或 2
            c0 = R / K;              % 中心占比
            c1 = 1 - c0;             % 圆环占比
            % 目标极角 theta：cos(theta) = (β - c0)/c1，夹紧到 [-1,1]
            ct = (beta - c0) / max(c1,1e-12); ct = max(-1,min(1,ct));
            st = sqrt(max(0,1-ct^2));  % sin(theta)

            cols = {};
            % --- 环：Q 组方位对 ---
            for q = 1:Q
                psi = 2*pi*(q-1)/Q;      % 均匀方位
                w   = cos(psi)*e1 + sin(psi)*e2;  % 圆环切向单位矢
                vA  = ct*n + st*w;       vA = vA / norm(vA);
                vB  = ct*n - st*w;       vB = vB / norm(vB);
                pA  = stokes_to_jones_unit(vA);   % Jones（全局相位固定：第一分量实非负）
                pB  = stokes_to_jones_unit(vB);
                % 两相位：{0, π/2}，确保 Hbar=0
                cols{end+1} = pA;        %#ok<*AGROW>
                cols{end+1} = 1j*pA;
                cols{end+1} = pB;
                cols{end+1} = 1j*pB;
            end
            % --- 若 K≡2 (mod 4)，补“中心” n（两相位），不破坏横向相消 ---
            if R==2
                pC = stokes_to_jones_unit(n);
                cols{end+1} = pC;
                cols{end+1} = 1j*pC;
            end

            pk = [cols{:}];
            if size(pk,2) > polar_num, pk = pk(:,1:polar_num); end
            pk = norm_cols(pk);
            pk_rad = jones_to_gamma_eta(pk);

            % === 4) 自检（不 backproject）：G 是否等于 G*，κR 是否 ≤ δ ===
            Gemp = (pk*pk')/size(pk,2); Gemp=(Gemp+Gemp')/2; Gemp = Gemp / max(real(trace(Gemp)),1e-12);
            T    = S*Gemp*S';  Tn = ((T+T')/2) / max(real(trace(T)),1e-12);
            kR   = norm(Tn - 0.5*eye(2), 2);
            if kR > opts.deltaR + 5e-6
                warning('κR=%.4g 超出 δ=%.4g（检查 K 与 δ；或放宽 δ/增大 K）', kR, opts.deltaR);
            end

        case 12  % === track12：已知 S,pj 的定向最优 Gbar ===
            % 说明：
            %   - 要求 opts.pj 为 2x1 Jones 目标极化；
            %   - 可选 opts.mode='pure'|'phase'|'tilt'（默认 'pure'）
            %   - 可选 opts.theta_deg（tilt 小角度，默认 6 度）
            %   - 可选 opts.use_pinv=true 时用 pinv 代替反演
            %   - 本 case 为“最优定向”设计，不再做 S-反投影（backproject 被忽略）
            if isempty(S)
                error('case 12 需要提供 S (2x2)。');
            end
            if ~isfield(opts,'pj') || isempty(opts.pj)
                error('case 12 需要在 opts.pj 中提供 2x1 的目标 Jones 向量。');
            end
            pj_local = opts.pj;
            if ~isequal(size(pj_local), [2,1])
                error('opts.pj 必须是 2x1 列向量。');
            end

            % 读取可选参数
            if ~isfield(opts,'mode'),       opts.mode = 'pure'; end
            if ~isfield(opts,'theta_deg'),  opts.theta_deg = 6; end
            if ~isfield(opts,'use_pinv'),   opts.use_pinv = false; end
            Mode = lower(opts.mode);
            theta = opts.theta_deg * pi/180;
            K = polar_num;

            % q = S^{-1} p_j 或 pinv
            if ~opts.use_pinv
                q = S \ pj_local;            % 数值更稳的左除
                if ~all(isfinite(q))
                    q = pinv(S) * pj_local;
                end
            else
                q = pinv(S) * pj_local;
            end
            nq = norm(q);
            if nq < 1e-12
                error('S^{-1} p_j 近似为零，S 可能奇异或 pj 异常。');
            end
            q0 = q / nq;                     % 单位化的最优方向

            % 生成轨道 pk
            pk = zeros(2, K);
            switch Mode
                case 'pure'
                    % 全部同向（或同向+任意全局相位）；Gbar = q0 q0'
                    for k = 1:K
                        pk(:,k) = q0;
                    end

                case 'phase'
                    % 相位匀旋；Gbar 同样等于 q0 q0'
                    phi = 2*pi*(0:K-1)/K;
                    for k = 1:K
                        pk(:,k) = q0 * exp(1j*phi(k));
                    end

                case 'tilt'
                    % 小角度向正交方向抖动：Gbar 为两纯态凸组合，增强稳健性
                    t0 = [-conj(q0(2)); conj(q0(1))];  % 与 q0 正交的单位向量
                    t0 = t0 / max(norm(t0),1e-12);
                    phi = 2*pi*(0:K-1)/K;
                    c = cos(theta); s = sin(theta);
                    for k = 1:K
                        pk(:,k) = c*q0 + s*exp(1j*phi(k))*t0;
                        pk(:,k) = pk(:,k) / max(norm(pk(:,k)),1e-12);
                    end

                otherwise
                    error('opts.mode="%s" 不支持（可选：pure|phase|tilt）。', Mode);
            end

            % 最终单位化，并输出 (gamma, eta)
            pk = norm_cols(pk);
            % 注意：本 case 已是 Tx 最优方向，不再做 backproject（会破坏最优性）
            pk_rad = jones_to_gamma_eta(pk);
 
            
        case 13  % OPT-RxRelax-v2: 约束 kappaR ≤ δ，最小化 kappaT（搜 λ 和方向 u）
            if isempty(S), error('case 13 需要 S'); end
            if ~isfield(opts,'deltaR'),     opts.deltaR = 0.03; end    % δ∈[0,0.03]
            if ~isfield(opts,'oriN'),       opts.oriN   = max(600, 120*polar_num); end  % 方向候选数
            if ~isfield(opts,'rng_seed'),   opts.rng_seed = 1; end
            
            K = polar_num;
            if mod(K,2)~=0
                warning('K 需为偶数以成对实现；已改为 K-1'); K = K-1;
            end
            
            % === 方向候选（覆盖整个球面） ===
            Uxyz = fibonacci_sphere(opts.oriN);      % 3xN
            [Uj, ~] = stokes3_to_jones(Uxyz);       % 2xN 纯态 Jones 方向
            Nori = size(Uj,2);
            
            best = struct('ok',false,'n2',[],'lam',[],'u',[],'kR',[],'kT',[],'T',[]);
            for n2 = 0:2:K
                lam = n2 / K;   % λ = 份额（trace(T)=1）
                kR_min = inf;   u_best = []; T_best = []; kT_best = inf;
                
                for i = 1:Nori
                    u = Uj(:,i); u = u / max(norm(u),1e-12);
                    up = [-conj(u(2)); conj(u(1))];   % 正交
                    Ttry = lam*(u*u') + (1-lam)*(up*up');       % Hermitian, trace=1
                    
                    Rtry = S*Ttry*S';
                    Rh   = (Rtry+Rtry')/2;
                    Rnh  = Rh / max(real(trace(Rh)),1e-12);
                    kR   = norm(Rnh - 0.5*eye(2), 2);
                    if kR < kR_min
                        kR_min = kR; u_best = u; T_best = Ttry;
                    end
                end
                
                if kR_min <= opts.deltaR + 1e-12
                    % 该 λ 可行：其 kappaT 仅由 λ 决定
                    kT = abs(lam - 0.5);
                    if ~best.ok || (kT < best.kT - 1e-12) || (abs(kT-best.kT)<=1e-12 && kR_min < best.kR)
                        best.ok = true; best.n2 = n2; best.lam = lam; best.u = u_best;
                        best.kR = kR_min; best.kT = kT; best.T = T_best;
                    end
                end
            end
            
            if ~best.ok
                error('没有满足 kappaR ≤ δ=%.3f 的可行解（K=%d 太小或 δ 过小）。建议增大 K 或放宽 δ。', opts.deltaR, polar_num);
            end
            
            % === 用 {0,π} 成对相位精确实现 T(best) ===
            n2 = best.n2; n1 = K - n2;
            u  = best.u;  up = [-conj(u(2)); conj(u(1))];
            pk = zeros(2,K); ptr=0;
            
            % λ 对应方向 u 放 n2 段（成对）
            for i=1:(n2/2)
                ptr=ptr+1; pk(:,ptr) = u;
                ptr=ptr+1; pk(:,ptr) = -u;
            end
            % 1-λ 对应方向 up 放 n1 段（成对）
            for i=1:(n1/2)
                ptr=ptr+1; pk(:,ptr) = up;
                ptr=ptr+1; pk(:,ptr) = -up;
            end
            
            % 二次相位分组：一半对整体加 π/2，使 M2≈0（不影响 T/R）
            half_pairs = K/4;
            for j=1:half_pairs
                idx1 = 2*j-1; idx2 = 2*j;           % 一对
                pk(:,idx1) = pk(:,idx1)*1j;         % 乘 e^{jπ/2}
                pk(:,idx2) = pk(:,idx2)*1j;
            end
            
            pk = norm_cols(pk);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 14

            M = K;
            p_hat = opts.pj;
            betas.beta1 = 0.10;        % 目标沿干扰方向能量泄露上限
            betas.beta2 = 0.05;        % 干扰自耦合系数上限
            betas.beta3 = 0.05;        % 交叉项实部上限
            opts = struct();
            opts.method = 'grid';      % 无CVX就用默认网格-细化
            opts.alpha  = 1;           % JSR 归一化
            opts.sigma2 = 1;           % 噪声归一化
            
            [pk, ~, ~, ~] = build_track14(S, p_hat, M, betas, opts);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 15 % H单极化
            gamma = ones(1, K) .* (pi / 2);
            eta   = zeros(1, K);
            pk    = [cos(gamma); sin(gamma) .* exp(1j * eta)];
            pk    = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
        case 16 % 圆极化
            gamma = ones(1, K) .* pi / 4;
            eta   = ones(1, K) .* pi / 2;
            gamma = zeros(1, K);
            eta   = zeros(1, K);
            pk    = [cos(gamma); sin(gamma) .* exp(1j * eta)];
            pk    = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        otherwise
            error('未知 polar_type=%d', polar_type);
    end
end

% ===================== 工具函数 =====================

function X = norm_cols(X)
    n = sqrt(sum(abs(X).^2,1)); n(n<1e-12)=1; X = bsxfun(@rdivide,X,n);
end

function pk = stokes3_to_jones_and_backproj(U, S, do_backproj)
    % U: 3xK, (x,y,z) on unit sphere; RC 在 +Z
    x = U(1,:); y = U(2,:); z = U(3,:);
    gamma = 0.5*acos(max(-1,min(1,x)));
    eta   = atan2(z, y);  % RC +Z 约定
    uk    = [cos(gamma); sin(gamma).*exp(1j*eta)];
    uk    = norm_cols(uk);
    if do_backproj && ~isempty(S)
        % S^? = S^H (S S^H)^{-1}
        pk = S' * ((S*S') \ uk);
        pk = norm_cols(pk);
    else
        pk = uk;
    end
end

function pk_rad = jones_to_gamma_eta(pk)
    pk = norm_cols(pk);
    a  = pk(1,:); b = pk(2,:);
    gamma = atan2(abs(b), abs(a));
    eta   = angle(b) - angle(a);
    eta   = wrap_to_pi_local(eta);
    pk_rad = [gamma; eta];
end

function U = fibonacci_sphere(K)
    golden = (1 + sqrt(5))/2;
    i = (0:K-1) + 0.5;
    z = 1 - 2*i./K;
    r = sqrt(max(0,1 - z.^2));
    theta = 2*pi*i/golden;
    x = r.*cos(theta); 
    y = r.*sin(theta);
    U = [x; y; z];
    U = norm_cols(U);
end

function U = icosahedron_vertices()
    % 返回大约 20 个均匀点（12 顶点 + 若干边中点），列单位向量
    phi = (1+sqrt(5))/2;
    V = [ ...
         0  1  phi;
         0 -1  phi;
         0  1 -phi;
         0 -1 -phi;
         1  phi  0;
        -1  phi  0;
         1 -phi  0;
        -1 -phi  0;
         phi  0  1;
        -phi  0  1;
         phi  0 -1;
        -phi  0 -1]';
    V = norm_cols(V);
    % 简单补点使分布更圆
    E = [ (V(:,1)+V(:,5))/2, (V(:,2)+V(:,5))/2, (V(:,1)+V(:,9))/2, ...
          (V(:,3)+V(:,11))/2, (V(:,4)+V(:,12))/2, (V(:,6)+V(:,10))/2, ...
          (V(:,7)+V(:,11))/2, (V(:,8)+V(:,12))/2 ];
    E = norm_cols(E);
    U = [V, E];
    U = norm_cols(U);
end

function x = wrap_to_pi_local(x)
    % 无 Mapping Toolbox 时自带的 wrapToPi
    x = mod(x + pi, 2*pi) - pi;
end


function [J, rad] = stokes3_to_jones(U)
    x=U(1,:); y=U(2,:); z=U(3,:);
    gamma = 0.5*acos(max(-1,min(1,x)));
    eta   = atan2(z,y);
    J = [cos(gamma); sin(gamma).*exp(1j*eta)];
    J = norm_cols(J);
    if nargout>1, rad=[gamma;eta]; end
end

function [u_best, lam_best] = solve_u_lambda_minKappaT(S, delta, oriN)
% 网格搜索 u（Jones 纯态方向）+ 离散 λ = n2/K 近似，
% 返回满足 κR≤δ 时使 κT 最小的 (u,λ)。若多解并列，取 κR 最小者。
    Uxyz = fibonacci_sphere(max(600, oriN));        % 3xN 候选方向（Stokes）
    [Uj, ~] = stokes3_to_jones(Uxyz);               % 2xN 纯 Jones
    lam_grid = linspace(0, 0.5, 201);               % 细网格（可根据 K 调整）
    best = struct('ok',false,'u',[],'lam',[],'kT',inf,'kR',inf);
    for i = 1:size(Uj,2)
        u = Uj(:,i) / max(norm(Uj(:,i)),1e-12);
        up = [-conj(u(2)); conj(u(1))];
        for lam = lam_grid
            G = (1-lam)*(u*u') + lam*(up*up'); G=(G+G')/2; G=G/max(real(trace(G)),1e-12);
            T = S*G*S'; Tn=((T+T')/2)/max(real(trace(T)),1e-12);
            kR = norm(Tn - 0.5*eye(2), 2);
            if kR <= delta + 1e-12
                kT = abs(lam - 0.5);
                if (~best.ok) || (kT < best.kT-1e-12) || (abs(kT-best.kT)<=1e-12 && kR < best.kR)
                    best.ok=true; best.u=u; best.lam=lam; best.kT=kT; best.kR=kR;
                end
            end
        end
    end
    if ~best.ok
        error('δ=%.3f 下无可行解（可增大 K 或放宽 δ）。', delta);
    end
    u_best = best.u; lam_best = best.lam;
end

function p = stokes_to_jones_unit(v)
% v∈S^2 -> 规范化 Jones（第一分量取非负实数的相位规范）
    x=v(1); y=v(2); z=v(3);
    gamma = 0.5*acos(max(-1,min(1,x)));
    eta   = atan2(z, y);
    p     = [cos(gamma); sin(gamma).*exp(1j*eta)];
    % 统一全局相位：使第一分量实且非负
    ph = angle(p(1)); p = p*exp(-1j*ph);
    if real(p(1))<0, p = -p; end
end



