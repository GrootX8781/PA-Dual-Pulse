function info = pk_info(pk, S, varargin)
    
    pk = norm_cols(pk);
    K  = size(pk, 2);
    
    Tbar = zeros(2, 2);
    for k = 1:K
        Tbar = Tbar + pk(:, k) * pk(:, k)';
    end
    Tbar = Tbar / K;
    
    rk = S * pk;
%     rk = norm_cols(rk);
    Rbar = zeros(2, 2);
    e    = zeros(2, k);
    for k = 1:K
        Rbar = Rbar + rk(:, k) * rk(:, k)';
        e(k) = real(rk(:, k)' * rk(:, k));
    end
    Rbar = Rbar / K;
    
    % ∑¢…‰
    Th       = (Tbar + Tbar') / 2;
    trT      = max(real(trace(Th)), 1e-12);
    Th       = Th / trT;
    lamT     = sort(real(eig(Th)), 'ascend');
    lammin_T = max(lamT(1), 0);
    lammax_T = max(lamT(2), 0);
    kappa_T  = norm(Th - 0.5 * eye(2), 2);

    phi = zeros(1, K);
    for k = 1:K
        a1 = pk(1, k);
        b1 = pk(2, k);
        if abs(a1) >= 1e-8
            phi(k) = angle(a1);
        else
            phi(k) = angle(b1);
        end
    end
    M1 = abs(sum(exp(1j * phi))) / K;
    M2 = abs(sum(exp(1j * 2 * phi))) / K;
    
    % Ω” ’
    Rh       = (Rbar + Rbar') / 2;
    trR      = max(real(trace(Rh)), 1e-12);
    Rh       = Rh / trR;
    lamR     = sort(real(eig(Rh)), 'ascend');
    lammin_R = max(lamR(1), 0);
    lammax_R = max(lamR(2), 0);
    kappa_R  = norm(Rh - 0.5 * eye(2), 2);
         
    info = struct();
    info.T       = Tbar;
    info.lamT    = lamT;
    info.lamminT = lammin_T;
    info.lammaxT = lammax_T;
    info.kappaT  = kappa_T;
    info.M1      = M1;
    info.M2      = M2;
    
    info.R       = Rbar;
    info.lamR    = lamR;
    info.lamminR = lammin_R;
    info.lammaxR = lammax_R;
    info.kappaR  = kappa_R;

end

function X = norm_cols(X)
    n = sqrt(sum(abs(X).^2, 1)); 
    n(n < 1e-12) = 1; 
    X = bsxfun(@rdivide, X, n);
end

