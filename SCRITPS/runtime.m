clc; clear; close all;
rng(1);

% -----------------------------
% Basic setup
% -----------------------------
cvx_solver sdpt3;
cvx_quiet(true);

n_warmup = 3;
n_repeat = 50;

delta_max = 0.20;
beta1 = 0.10;
beta2 = 0.05;
beta3 = 0.05;

% Representative PSMs
S_list = cell(3,1);
S_list{1} = S1;
S_list{2} = S2;
S_list{3} = S3;

S_name = {'S1','S2','S3'};

% One nominal estimated jamming polarization
phat_gamma = 30*pi/180;
phat_eta   = 40*pi/180;
phat = [cos(phat_gamma);
        sin(phat_gamma)*exp(1j*phat_eta)];
phat = phat / norm(phat);

fprintf('===============================================\n');
fprintf('Benchmark of eqopt1 and eqopt2 using CVX\n');
fprintf('Solver: SDPT3\n');
fprintf('Warm-up runs: %d\n', n_warmup);
fprintf('Measured runs: %d\n', n_repeat);
fprintf('===============================================\n\n');

for is = 1:numel(S_list)
    S = S_list{is};
    fprintf('---- %s ----\n', S_name{is});

    % warm-up
    for k = 1:n_warmup
        solve_eqopt1_cvx(S, delta_max);
        solve_eqopt2_cvx(S, phat, beta1, beta2, beta3);
    end

    % eqopt1 timing
    t1 = zeros(n_repeat,1);
    status1 = strings(n_repeat,1);
    for k = 1:n_repeat
        tic;
        [~, status1(k)] = solve_eqopt1_cvx(S, delta_max);
        t1(k) = toc * 1e3; % ms
    end

    % eqopt2 timing
    t2 = zeros(n_repeat,1);
    status2 = strings(n_repeat,1);
    for k = 1:n_repeat
        tic;
        [~, status2(k)] = solve_eqopt2_cvx(S, phat, beta1, beta2, beta3);
        t2(k) = toc * 1e3; % ms
    end

    ok1 = ismember(status1, ["Solved","Inaccurate/Solved"]);
    ok2 = ismember(status2, ["Solved","Inaccurate/Solved"]);

    fprintf('eqopt1 success: %d / %d\n', sum(ok1), n_repeat);
    fprintf('eqopt1 mean   : %.4f ms\n', mean(t1(ok1)));
    fprintf('eqopt1 std    : %.4f ms\n', std(t1(ok1)));

    fprintf('eqopt2 success: %d / %d\n', sum(ok2), n_repeat);
    fprintf('eqopt2 mean   : %.4f ms\n', mean(t2(ok2)));
    fprintf('eqopt2 std    : %.4f ms\n\n', std(t2(ok2)));
end

% ==========================================
% Local functions
% ==========================================
function [G, status] = solve_eqopt1_cvx(S, delta_max)
    cvx_begin sdp quiet
        variable G(2,2) hermitian
        G == hermitian_semidefinite(2);
        minimize( norm(G - 0.5*eye(2), 2) )
        subject to
            trace(G) == 1;
            norm(S*G*S' - 0.5*eye(2), 2) <= delta_max;
    cvx_end
    status = string(cvx_status);
end

function [G, status] = solve_eqopt2_cvx(S, phat, beta1, beta2, beta3)
    Pperp = eye(2) - phat*phat';
    C = S' * Pperp * S;
    B = S' * phat * phat' * S;
    D = conj(phat) * phat.';       % overline(phat) * phat^T
    M = conj(phat) * phat' * S;    % overline(phat) * phat^H * S

    cvx_begin sdp quiet
        variable G(2,2) hermitian
        G == hermitian_semidefinite(2);
        maximize( real(trace(C*G)) )
        subject to
            trace(G) == 1;
            real(trace(B*G)) <= beta1;
            real(trace(D*G)) <= beta2;
            abs(real(trace(M*G))) <= beta3;
    cvx_end
    status = string(cvx_status);
end