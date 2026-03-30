
function v = get_max_ev(Y, lambda_scale)
    if nargin < 2
        R = (Y * Y') / max(1, size(Y, 2)); 
        R = (R + R')/2; 
        [V, D] = eig(R);
    else
        R = (Y * Y') / max(1, size(Y, 2));
        R = (R + R')/2;
        lambda = lambda_scale * trace(R) / 2;
        Rdl = R + lambda * eye(2);
        [V, D] = eig((Rdl + Rdl') / 2);
    end
    [~, ix] = max(real(diag(D)));
    v = V(:, ix);
    v = v / max(norm(v), eps);
end