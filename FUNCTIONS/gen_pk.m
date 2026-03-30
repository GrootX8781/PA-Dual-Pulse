function [pk, pk_rad] = gen_pk(polar_num, type, S, opts)

    if nargin < 3 || isempty(S)
        S = [];
    end
    if nargin < 4
        opts = struct();
    end
    if ~isfield(opts, 'backproject')
        opts.backproject = true;
    end
    
    M = polar_num;
    do_backproj = opts.backproject;
    
    switch type
        case 1 % 随机
            g  = sort(rand(1, M) * (pi / 2));
            e  = (rand(1, M) - 0.5) * (2 * pi);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
        case 2 % 变极化 H极化-V极化
            g  = linspace(0, pi / 2, M);
            e  = zeros(1, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 3 % 变极化 RC极化-LC极化
            g  = (pi / 4) * ones(1, M);
            e  = linspace(pi / 2, -pi / 2, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 4 % 变极化 γη均匀扫描
            g  = linspace(0, pi/2, M);
            e  = linspace(-pi, pi, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);

        case 5 % EQR 子午环
            g  = (pi / 4) * ones(1, M);
            e  = linspace(0, 2 * pi, M + 1); e(end) = [];
            U  = [cos(2 * g); sin(2 * g) .* cos(e); sin(2 * g) .* sin(e)];
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 6 % TRG 双环
            if ~isfield(opts, 'trg_delta')
                opts.trg_delta = 15;
            end
            d    = opts.trg_delta * pi / 180;
            phi  = 0.7 * (pi / (M - floor(M / 2)));
            M1   = floor(M / 2);
            M2   = M - M1;
            e1 = linspace(0, 2 * pi, M1 + 1); e1(end) = [];
            e2 = linspace(phi, 2 * pi + phi, M2 + 1); e2(end) = [];
            g1 = (pi / 4 - d) * ones(1, M1);
            g2 = (pi / 4 + d) * ones(1, M2);
            U1 = [cos(2 * g1); sin(2 * g1) .* cos(e1); sin(2 * g1) .* sin(e1)];
            U2 = [cos(2 * g2); sin(2 * g2) .* cos(e2); sin(2 * g2) .* sin(e2)];
            U  = [U1, U2];
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 7 % ICO 正二十面体
            U0 = icosahedron_vertices();
            n  = size(U0, 2);
            if M <= n
                idx = round(linspace(1, n, M));
                idx(idx < 1) = 1; idx(idx > n) = n;
                idx = unique(idx, 'stable');
                if numel(idx) < M
                    pool = setdiff(1:n, idx, 'stable');
                    need = M - numel(idx);
                    idx  = [idx, pool(1:need)];
                end
                U = U0(:, idx);
            else
                rep = ceil(M / n);
                U = repmat(U0, 1, rep);
                U = U(:, 1 : M);
            end
            U  = norm_cols(U);
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 8 % FIB fibonacci
            U = fibonacci_sphere(M);
            pk = stokes3_to_jones_and_backproj(U, S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 9 % COV
            load('best_cover.mat');
            if M < 4 || M > 30
                error('WRONG POLAR NUM');
            end
            U = cell2mat(best_cover(M - 3));
            pk = stokes3_to_jones_and_backproj(U.', S, opts.backproject);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 10 % OPT1 干扰极化估计轨道
  
            pk = build_opt1(S, M, opts);
%             pk = build_opt1_plus(S, M, opts);
           
            pk = norm_cols(pk);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 11 % OPT2 bestpk best sinr轨道
            phat = opts.opt2_pj;
            betas.beta1 = 0.10;
            betas.beta2 = 0.05;
            betas.beta3 = 0.05;
            
            pk = build_opt2(S, phat, M, betas);
            pk_rad = jones_to_gamma_eta(pk);
            
        case 12 % H单极化
            g  = zeros(1, M);
            e  = zeros(1, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
        case 13 % V单极化
            g  = ones(1, M) .* (pi / 2);
            e  = zeros(1, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
        case 14
            g = ones(1, M) .* (pi / 4);
            e = zeros(1, M);
            pk = [cos(g); sin(g) .* exp(1j * e)];
            pk = norm_cols(pk);
            if do_backproj && ~isempty(S)
                pk = S' * ((S * S') \ pk);
                pk = norm_cols(pk);
            end
            pk_rad = jones_to_gamma_eta(pk);
            
    end
    
end

%%
function X = norm_cols(X)
    n = sqrt(sum(abs(X) .^ 2, 1)); 
    n(n < 1e-12) = 1; 
    X = bsxfun(@rdivide, X, n);
end

function pk_rad = jones_to_gamma_eta(pk)
    pk = norm_cols(pk);
    a  = pk(1, :); 
    b  = pk(2, :);
    gamma = atan2(abs(b), abs(a));
    eta   = angle(b) - angle(a);
    eta   = wrapToPi(eta);
    pk_rad = [gamma; eta];
end

function [J, rad] = stokes3_to_jones(U)
    x = U(1, :); y = U(2, :); z = U(3, :);
    gamma = 0.5 * acos(max(-1, min(1, x)));
    eta   = atan2(z, y);
    J = [cos(gamma); sin(gamma) .* exp(1j * eta)];
    J = norm_cols(J);
    if nargout>1, rad=[gamma; eta]; end
end

function pk = stokes3_to_jones_and_backproj(U, S, do_backproj)
    x = U(1, :); 
    y = U(2, :); 
    z = U(3, :);
    gamma = 0.5 * acos(max(-1, min(1, x)));
    eta   = atan2(z, y);
    uk    = [cos(gamma); sin(gamma) .* exp(1j * eta)];
    uk    = norm_cols(uk);
    if do_backproj && ~isempty(S)
        pk = S' * ((S * S') \ uk);
        pk = norm_cols(pk);
    else
        pk = uk;
    end
end

function U = icosahedron_vertices()
    phi = (1 + sqrt(5)) / 2;
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
    E = [ (V(:,1)+V(:,5))/2, (V(:,2)+V(:,5))/2, (V(:,1)+V(:,9))/2, ...
          (V(:,3)+V(:,11))/2, (V(:,4)+V(:,12))/2, (V(:,6)+V(:,10))/2, ...
          (V(:,7)+V(:,11))/2, (V(:,8)+V(:,12))/2 ];
    E = norm_cols(E);
    U = [V, E];
    U = norm_cols(U);
end

function U = fibonacci_sphere(M)
    golden = (1 + sqrt(5)) / 2;
    i = (0:M - 1) + 0.5;
    z = 1 - 2 * i ./ M;
    r = sqrt(max(0, 1 - z.^2));
    theta = 2 * pi * i / golden;
    x = r .* cos(theta); 
    y = r .* sin(theta);
    U = [x; y; z];
    U = norm_cols(U);
end

function pk = build_opt1(S, M, opts)
    if isempty(S), error('case 10 需要 S'); end
    if ~isfield(opts, 'opt1_deltaR'), opts.opt1_deltaR = 0.03; end

    if mod(M, 2) ~= 0
        warning('M needs to be odd, switch to M - 1');
        M = M - 1;
    end

    oriN = max(600, 120 * M);
    U = fibonacci_sphere(oriN);
    [Uj, ~] = stokes3_to_jones(U);
    Nori = size(Uj, 2);

    best = struct('ok', false, 'n2', [], 'lam', [], ...
        'u', [], 'kR', [], 'kT', [], 'T', []);

    for n2 = 0:2:M
        lam = n2 / M;
        kR_min = inf; u_best = []; T_best = [];
        
        for i = 1:Nori
            u  = Uj(:, i); u = u / max(norm(u), 1e-12);
            up = [-conj(u(2)); conj(u(1))];
            Ttry = lam * (u * u') + (1 - lam) * (up * up');
            
            Rtry = S * Ttry * S';
            Rh   = (Rtry + Rtry') / 2;
            Rnh  = Rh / max(real(trace(Rh)), 1e-12);
            kR   = norm(Rnh - 0.5*eye(2), 2);
            if kR < kR_min
                kR_min = kR; u_best = u; T_best = Ttry;
            end
        end
        
        if kR_min <= opts.opt1_deltaR + 1e-12
            kT = abs(lam - 0.5);
            if ~best.ok || (kT < best.kT - 1e-12) || (abs(kT-best.kT) <= 1e-12 && kR_min < best.kR)
                best.ok = true; best.n2 = n2; best.lam = lam; best.u = u_best;
                best.kR = kR_min; best.kT = kT; best.T = T_best;
            end
        end
    end
    
    if ~best.ok
        error('NO OPTIMAL TRACK UNDER δ=%.3f\n', opts.opt1_deltaR);
    end
    
    n2 = best.n2;     n1  = M - n2;
    u  = best.u;      up  = [-conj(u(2)); conj(u(1))];
    pk = zeros(2, M); ptr = 0;
    
    for i = 1:(n2 / 2)
        ptr = ptr + 1; pk(:, ptr) = u;
        ptr = ptr + 1; pk(:, ptr) = -u;
    end
    for i = 1:(n1 / 2)
        ptr = ptr + 1; pk(:, ptr) = up;
        ptr = ptr + 1; pk(:, ptr) = -up;
    end
    half_pairs = M / 4;
    for j = 1:half_pairs
        idx1 = 2 * j - 1; idx2 = 2 * j;
        pk(:, idx1) = pk(:, idx1) * 1j;
        pk(:, idx2) = pk(:, idx2) * 1j;
    end

end


function pk = build_opt2(S, phat, M, betas)

    beta1 = betas.beta1; 
    beta2 = betas.beta2; 
    beta3 = betas.beta3;

    polar_num = M;
    
    Pperp = eye(2) - phat * phat';
    C     = S' * Pperp * S;
    B     = S' * phat * (phat') * S;
    D     = conj(phat) * phat';
    M     = conj(phat) * (phat') * S;
    
    best = -inf; pbest = [1; 0];
    Ng = 181;
    Ne = 360;
    gamma = linspace(0, pi / 2, Ng);
    eta   = linspace(0, 2 * pi, Ne);
    for g = gamma
        cg = cos(g); sg = sin(g);
        for e = eta
            p = [cg; sg * exp(1j * e)];
            if real(p' * B * p) > beta1, continue; end
            if real(p' * D * p) > beta2, continue; end
            mu = real(p' * M * p);
            if abs(mu) > beta3, continue; end
            val = real(p' * C * p);
            if val > best, best = val; pbest = p; end
        end
    end
    
    for it = 1:3
        g0 = acos(abs(pbest(1))); e0 = angle(pbest(2));
        dg = (pi / 2) / max(30 * 2 ^ it, 60);
        de = (2 * pi) / max(60 * 2 ^ it, 120);
        gg = g0 + (-6:6) * dg;
        ee = e0 + (-12:12) * de;
        for g = gg
            if g < 0 || g > pi / 2, continue; end
            cg = cos(g); sg = sin(g);
            for e = ee
                p = [cg; sg * exp(1j * e)];
                if real(p' * B * p) > beta1, continue; end
                if real(p' * D * p) > beta2, continue; end
                mu = real(p' * M * p);
                if abs(mu) > beta3, continue; end
                val = real(p' * C * p);
                if val > best, best = val; pbest = p; end
            end
        end
        pbest = pbest / max(norm(pbest),1e-12);
        Gstar = pbest * pbest';
    end
    
    [V, Lam]  = eig((Gstar + Gstar') / 2);
    lam       = real(diag(Lam)); 
    [lam, ix] = sort(lam, 'descend'); 
    V = V(:, ix);
    
    if lam(1) >= 1-1e-6
        p1 = V(:, 1) * exp(-1j * angle(V(1, 1)));
        pk = repmat(p1, 1, polar_num);
    else
        p1 = V(:, 1) * exp(-1j * angle(V(1, 1)));
        p2 = V(:, 2) * exp(-1j * angle(V(1, 1)));
        m1 = max(round(lam(1) * polar_num), 0); 
        m1 = min(m1, polar_num); 
        m2 = polar_num - m1;
        pk = [repmat(p1, 1, m1), repmat(p2, 1, m2)];
    end

end
