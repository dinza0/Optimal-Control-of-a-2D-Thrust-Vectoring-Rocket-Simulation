%% Part 5.2 Controllability and Gramian Comparison
clear; clc; close all;

%% Parameters
p.g     = 9.81;      % m/s^2
p.I     = 10;        % kg*m^2
p.m0    = 20;        % kg
p.alpha = 0.005;     % kg/(N*s)

tf = 5.0;            % finite horizon, seconds
n  = 7;

%% Hover linearization, frozen mass
A_hover = makeA(p.g, p.m0);
B_hover = makeB(p.m0, p.I, p.alpha);

C_hover = controllabilityMatrix(A_hover, B_hover);
sC_hover = svd(C_hover);

fprintf('\n===== HOVER LINEARIZATION =====\n');
fprintf('rank(C_hover)     = %d\n', rank(C_hover));
fprintf('sigma_min(C_hover)= %.6e\n', min(sC_hover));
fprintf('sigma_max(C_hover)= %.6e\n', max(sC_hover));
fprintf('cond(C_hover)     = %.6e\n', max(sC_hover)/min(sC_hover));

[W_hover, eigW_hover] = gramianLTV(@(t) A_hover, @(t) B_hover, tf, n);

fprintf('\nHover Gramian:\n');
fprintf('lambda_min(W_hover)= %.6e\n', min(eigW_hover));
fprintf('lambda_max(W_hover)= %.6e\n', max(eigW_hover));
fprintf('cond(W_hover)      = %.6e\n', max(eigW_hover)/min(eigW_hover));

%% Weakly controllable time-varying descent linearization
% Nominal powered descent:
% theta_star approx 0, T_star = m_star*g
% m_dot_star = -alpha*m_star*g

mstar = @(t) p.m0 * exp(-p.alpha*p.g*t);

A_weak = @(t) makeA(p.g, mstar(t));
B_weak = @(t) makeB(mstar(t), p.I, p.alpha);

fprintf('\n===== WEAK TIME-VARYING DESCENT LINEARIZATION =====\n');

% Snapshot controllability matrices at several times
time_samples = [0, tf/2, tf];

for k = 1:length(time_samples)
    tk = time_samples(k);
    Ak = A_weak(tk);
    Bk = B_weak(tk);
    Ck = controllabilityMatrix(Ak, Bk);
    sCk = svd(Ck);

    fprintf('\nSnapshot t = %.2f s, m*(t)=%.4f kg\n', tk, mstar(tk));
    fprintf('rank(C_snapshot)      = %d\n', rank(Ck));
    fprintf('sigma_min(C_snapshot) = %.6e\n', min(sCk));
    fprintf('sigma_max(C_snapshot) = %.6e\n', max(sCk));
    fprintf('cond(C_snapshot)      = %.6e\n', max(sCk)/min(sCk));
end

% Finite-horizon LTV Gramian
[W_weak, eigW_weak] = gramianLTV(A_weak, B_weak, tf, n);

fprintf('\nWeak LTV Gramian:\n');
fprintf('lambda_min(W_weak)= %.6e\n', min(eigW_weak));
fprintf('lambda_max(W_weak)= %.6e\n', max(eigW_weak));
fprintf('cond(W_weak)      = %.6e\n', max(eigW_weak)/min(eigW_weak));

%% Horizontal subsystem Gramian
% This isolates the weak horizontal chain:
% tau -> omega -> theta -> vx -> px

A_horiz = [0 1 0 0;
           0 0 p.g 0;
           0 0 0 1;
           0 0 0 0];

B_horiz = [0;
           0;
           0;
           1/p.I];

C_horiz = controllabilityMatrix(A_horiz, B_horiz);
sC_horiz = svd(C_horiz);

[W_horiz, eigW_horiz] = gramianLTV(@(t) A_horiz, @(t) B_horiz, tf, 4);

fprintf('\n===== HORIZONTAL SUBSYSTEM =====\n');
fprintf('rank(C_horiz)      = %d\n', rank(C_horiz));
fprintf('sigma_min(C_horiz) = %.6e\n', min(sC_horiz));
fprintf('cond(C_horiz)      = %.6e\n', max(sC_horiz)/min(sC_horiz));

fprintf('\nHorizontal Gramian:\n');
fprintf('lambda_min(W_horiz)= %.6e\n', min(eigW_horiz));
fprintf('lambda_max(W_horiz)= %.6e\n', max(eigW_horiz));
fprintf('cond(W_horiz)      = %.6e\n', max(eigW_horiz)/min(eigW_horiz));

%% Optional plotting of Gramian eigenvalues
figure;
semilogy(eigW_hover, 'o-', 'LineWidth', 1.5); hold on;
semilogy(eigW_weak, 's-', 'LineWidth', 1.5);
grid on;
xlabel('Eigenvalue index');
ylabel('Gramian eigenvalue');
legend('Hover', 'Weak LTV descent');
title('Finite-Horizon Controllability Gramian Eigenvalues');

%% ================= LOCAL FUNCTIONS =================

function A = makeA(g, m)
    A = [0 0 1 0 0 0 0;
         0 0 0 1 0 0 0;
         0 0 0 0 g 0 0;
         0 0 0 0 0 0 -g/m;
         0 0 0 0 0 1 0;
         0 0 0 0 0 0 0;
         0 0 0 0 0 0 0];
end

function B = makeB(m, I, alpha)
    B = [0 0;
         0 0;
         0 0;
         1/m 0;
         0 0;
         0 1/I;
         -alpha 0];
end

function C = controllabilityMatrix(A, B)
    n = size(A,1);
    C = [];
    for k = 0:n-1
        C = [C, A^k*B];
    end
end

function [W, eigW] = gramianLTV(Afun, Bfun, tf, n)
    W0 = zeros(n,n);
    w0 = W0(:);

    odeFun = @(t,w) gramianODE(t, w, Afun, Bfun, n);

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    [~, wSol] = ode45(odeFun, [0 tf], w0, opts);

    W = reshape(wSol(end,:).', n, n);
    W = 0.5*(W + W.');     % enforce symmetry numerically
    eigW = sort(eig(W));
end

function dwdt = gramianODE(t, w, Afun, Bfun, n)
    W = reshape(w, n, n);
    A = Afun(t);
    B = Bfun(t);

    dW = A*W + W*A.' + B*B.';
    dwdt = dW(:);
end