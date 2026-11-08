%% Part III Phase A: Minimum-Time Ascent using Indirect Shooting
clear; clc; close all;

%% Parameters
p.g = 9.81;
p.I = 10;
p.mh = 20;
p.alpha = 0.005;
 
p.Tmin = 0;
p.Tmax = 300;
p.taumin = -30;
p.taumax = 30;

p.dTmin = p.Tmin - p.mh*p.g;
p.dTmax = p.Tmax - p.mh*p.g;

p.A = makeA(p.g,p.mh);
p.B = makeB(p.mh,p.I,p.alpha);

%% Initial state
x0 = [0;0;0;0;0;0;0];

%% Terminal target, mass-free version
% Try vzsky = 0, +2, -2 in separate runs

%xsky6 = [20;25;0;0;0;0];
%xsky6 = [0;10;0;0;0;0];
xsky6 = [10;15;0;0;0;0];

%% Solver setup
%eps_list = [20, 10, 5, 2, 1, 0.5, 0.2, 0.1];
%eps_list = [20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01];
eps_list = [20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002];

% Initial guess: costates and final time
zguess = [
    -1;     % lambda_px
    -1;     % lambda_pz
    -0.5;   % lambda_vx
    -1;     % lambda_vz
     0.1;   % lambda_theta
     0.1;   % lambda_omega
     0.1;   % lambda_m
     6.0    % tf guess [s]
];

options = optimoptions('fsolve', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 5e4, ...
    'MaxIterations', 2e3, ...
    'FunctionTolerance', 1e-8, ...
    'StepTolerance', 1e-8);

best_z = [];
best_resnorm = inf;
best_eps = NaN;

for k = 1:length(eps_list)
    p.eps = eps_list(k);

    fprintf('\nSolving with eps = %.4f\n', p.eps);

    residualFun = @(z) phaseA_residual_massfree(z,x0,xsky6,p);

    lb = [-1e3*ones(7,1); 0.5];   % lower bounds: costates, tf >= 0.5 s
    ub = [ 1e3*ones(7,1); 30];    % upper bounds: costates, tf <= 30 s
    
    options = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxFunctionEvaluations', 5e4, ...
        'MaxIterations', 2e3, ...
        'FunctionTolerance', 1e-8, ...
        'StepTolerance', 1e-8);
    
    [zsol,resnorm,residual,exitflag,output] = lsqnonlin(residualFun,zguess,lb,ub,options);
    
    fprintf('exitflag = %d, residual norm = %.6e\n', exitflag, norm(residual));
    
    if exitflag <= 0 || norm(residual) > 1e-3
        warning('Solver did not converge for eps = %.4f. Do not trust this solution.', p.eps);
    else
        zguess = zsol;
    end
    
    if norm(residual) < best_resnorm
        best_resnorm = norm(residual);
        best_z = zsol;
        best_eps = p.eps;
    end

end

finalResidual = phaseA_residual_massfree(zsol,x0,xsky6,p);

if norm(finalResidual) > 1e-3 || zsol(8) <= 0.5
    error('No valid Phase A solution found. Residual too large or tf invalid.');
end

%% Simulate final solution

zsol = best_z;
p.eps = best_eps;

fprintf('\nUsing best epsilon = %.6f\n', p.eps);

lambda0 = zsol(1:7);
tf = zsol(8);

y0 = [x0; lambda0];

[t,y] = ode45(@(t,y) canonical_ode(t,y,p), [0 tf], y0);

x = y(:,1:7);
lam = y(:,8:14);

u = zeros(length(t),2);
H = zeros(length(t),1);
phiT = zeros(length(t),1);
phitau = zeros(length(t),1);

for i = 1:length(t)
    [ui,H(i),phiT(i),phitau(i)] = control_and_H(x(i,:).',lam(i,:).',p);
    u(i,:) = ui.';
end

Tphys = p.mh*p.g + u(:,1);

dT_bb = zeros(length(t),1);
tau_bb = zeros(length(t),1);

for i = 1:length(t)
    if phiT(i) < 0
        dT_bb(i) = p.dTmax;
    else
        dT_bb(i) = p.dTmin;
    end

    if phitau(i) < 0
        tau_bb(i) = p.taumax;
    else
        tau_bb(i) = p.taumin;
    end
end

T_bb = p.mh*p.g + dT_bb;

figure;
subplot(2,1,1);
plot(t,Tphys,'LineWidth',1.5); hold on;
plot(t,T_bb,'--','LineWidth',1.5);
grid on;
ylabel('T [N]');
legend('Smoothed','Bang-bang from switching function');
title('Smoothed vs Bang-Bang Thrust');

subplot(2,1,2);
plot(t,u(:,2),'LineWidth',1.5); hold on;
plot(t,tau_bb,'--','LineWidth',1.5);
grid on;
ylabel('\tau [N*m]');
xlabel('Time [s]');
legend('Smoothed','Bang-bang from switching function');
title('Smoothed vs Bang-Bang Torque');

fprintf('\nFinal time tf = %.6f s\n', tf);
fprintf('Terminal residual norm = %.6e\n', norm(phaseA_residual_massfree(zsol,x0,xsky6,p)));
fprintf('Hamiltonian final = %.6e\n', H(end));
fprintf('Hamiltonian initial = %.6e\n', H(1));
fprintf('Hamiltonian variation = %.6e\n', max(H)-min(H));
fprintf('Max physical thrust = %.6f N\n', max(Tphys));
fprintf('Min physical thrust = %.6f N\n', min(Tphys));
fprintf('Max abs torque = %.6f N*m\n', max(abs(u(:,2))));

%% Plots
figure;
plot(x(:,1),x(:,2),'LineWidth',1.8); grid on; axis equal;
xlabel('p_x [m]');
ylabel('p_z [m]');
title('Phase A Minimum-Time Ascent Trajectory');

figure;
subplot(2,1,1);
plot(t,Tphys,'LineWidth',1.5); grid on;
ylabel('T [N]');
title('Physical Thrust');

subplot(2,1,2);
plot(t,u(:,2),'LineWidth',1.5); grid on;
ylabel('\tau [N*m]');
xlabel('Time [s]');
title('Torque');

figure;
plot(t,H,'LineWidth',1.5); grid on;
xlabel('Time [s]');
ylabel('Hamiltonian');
title('Hamiltonian Along Phase A');

figure;
subplot(2,1,1);
plot(t,phiT,'LineWidth',1.5); grid on;
ylabel('\phi_T');
title('Thrust Switching Function');

subplot(2,1,2);
plot(t,phitau,'LineWidth',1.5); grid on;
ylabel('\phi_\tau');
xlabel('Time [s]');
title('Torque Switching Function');

%% ============================================================
%  Phase A: Exact Bang-Bang Switch-Time Parameterization
%  ============================================================

% This refines the smoothed tanh solution into an exact bang-bang solution.
% It uses the observed switching sequence:
% T:   Tmax -> Tmin
% tau: taumax -> taumin -> taumax -> taumin

fprintf('\n\n===== SWITCH-TIME PARAMETERIZATION =====\n');

% Initial guesses from smoothed solution / previous bang-bang result
tf0 = tf;

% Thrust switches: Tmax -> Tmin -> Tmax
tT10 = 1.90;
tT20 = 2.80;

% Torque switches: taumax -> taumin -> taumax -> taumin
ttau10 = 0.48;
ttau20 = 1.62;
ttau30 = 2.68;

% z = [lambda0(7); tf; tT1; tT2; ttau1; ttau2; ttau3]
zsw0 = [lambda0; tf0; tT10; tT20; ttau10; ttau20; ttau30];

lb = [-1e4*ones(7,1); 0.5; 0; 0; 0; 0; 0];
ub = [ 1e4*ones(7,1); 20; 20; 20; 20; 20; 20];

opts_sw = optimoptions('lsqnonlin', ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 1e5, ...
    'MaxIterations', 5e3, ...
    'FunctionTolerance', 1e-10, ...
    'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10);

resfun_sw = @(z) switch_time_residual(z,x0,xsky6,p);

[zsw,resnorm_sw,resid_sw,exitflag_sw,output_sw] = lsqnonlin( ...
    resfun_sw,zsw0,lb,ub,opts_sw);

fprintf('\nSwitch-time solve exitflag = %d\n', exitflag_sw);
fprintf('Switch-time residual norm = %.6e\n', norm(resid_sw));

% Simulate final bang-bang solution
[tbb,ybb,ubb,Hbb,phiTbb,phitaubb] = simulate_switch_solution(zsw,x0,p);

xbb = ybb(:,1:7);
lbb = ybb(:,8:14);

Tbb = p.mh*p.g + ubb(:,1);

fprintf('\n===== EXACT BANG-BANG RESULT =====\n');
fprintf('tf                       = %.6f s\n', zsw(8));
fprintf('T switches               = %.6f, %.6f s\n', zsw(9), zsw(10));
fprintf('tau switches             = %.6f, %.6f, %.6f s\n', zsw(11), zsw(12), zsw(13));
fprintf('terminal residual norm   = %.6e\n', norm(xbb(end,1:6).'-xsky6));
terminal_error = xbb(end,1:6).' - xsky6;

fprintf('\nTerminal error components:\n');
fprintf('px error     = %.6e\n', terminal_error(1));
fprintf('pz error     = %.6e\n', terminal_error(2));
fprintf('vx error     = %.6e\n', terminal_error(3));
fprintf('vz error     = %.6e\n', terminal_error(4));
fprintf('theta error  = %.6e\n', terminal_error(5));
fprintf('omega error  = %.6e\n', terminal_error(6));
fprintf('Hamiltonian final        = %.6e\n', Hbb(end));
fprintf('Hamiltonian initial      = %.6e\n', Hbb(1));
fprintf('Hamiltonian variation    = %.6e\n', max(Hbb)-min(Hbb));
fprintf('max physical thrust      = %.6f N\n', max(Tbb));
fprintf('min physical thrust      = %.6f N\n', min(Tbb));
fprintf('max abs torque           = %.6f N*m\n', max(abs(ubb(:,2))));

%% Plots for exact bang-bang solution

figure;
plot(xbb(:,1),xbb(:,2),'LineWidth',1.8);
grid on; axis equal;
xlabel('p_x [m]');
ylabel('p_z [m]');
title('Exact Bang-Bang Phase A Trajectory');

figure;
subplot(2,1,1);
plot(tbb,Tbb,'LineWidth',1.5);
grid on;
ylabel('T [N]');
title('Exact Bang-Bang Physical Thrust');

subplot(2,1,2);
plot(tbb,ubb(:,2),'LineWidth',1.5);
grid on;
ylabel('\tau [N*m]');
xlabel('Time [s]');
title('Exact Bang-Bang Torque');

figure;
plot(tbb,Hbb,'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('Hamiltonian');
title('Hamiltonian for Exact Bang-Bang Solution');

figure;
subplot(2,1,1);
plot(tbb,phiTbb,'LineWidth',1.5);
grid on;
ylabel('\phi_T');
title('Thrust Switching Function');

subplot(2,1,2);
plot(tbb,phitaubb,'LineWidth',1.5);
grid on;
ylabel('\phi_\tau');
xlabel('Time [s]');
title('Torque Switching Function');


%% ============================================================
%  Phase A Sensitivity Study: Terminal Vertical Velocity
%  ============================================================
% This block repeats the exact bang-bang switch-time solve for different
% terminal vertical velocity requirements at the sky rendezvous.
%
% IMPORTANT:
%   The sensitivity cases use the SAME assumed switching structure as the
%   final exact bang-bang solution above:
%       T:   Tmax -> Tmin -> Tmax
%       tau: taumax -> taumin -> taumax -> taumin
%
%   If a case does not converge to a small residual, that may mean the
%   terminal condition requires a different switching structure. In that
%   case, do not report that case as a valid exact solution without further
%   switch-structure refinement.

runTerminalVelocitySensitivity = true;

if runTerminalVelocitySensitivity

    fprintf('\n\n===== TERMINAL VELOCITY SENSITIVITY STUDY =====\n');

    % Change these values to test different terminal vertical velocities.
    % The baseline uses the original terminal velocity in xsky6.
    vz_sweep = [0; 2; -2];   % [m/s]

    case_names = { ...
        'Zero terminal velocity'; ...
        'Positive vertical terminal velocity'; ...
        'Negative vertical terminal velocity'};

    % Use a cell array instead of a struct array. This avoids MATLAB errors
    % caused by assigning structs that contain nested solver-output structs.
    sens = cell(length(vz_sweep),1);

    for j = 1:length(vz_sweep)
        xsky_case = xsky6;
        xsky_case(4) = vz_sweep(j);  % terminal v_z requirement

        if j == 1 && abs(vz_sweep(j) - xsky6(4)) < 1e-12
            % Reuse the already-computed exact bang-bang solution.
            sens{j} = package_switch_result(case_names{j},vz_sweep(j),zsw,exitflag_sw,resid_sw,xsky_case,x0,p);
        else
            % Use the baseline exact bang-bang solution as the initial guess.
            z0_case = zsw;

            % Small final-time adjustment for the initial guess only.
            % Positive terminal v_z usually needs less final braking; negative
            % terminal v_z may need more correction.
            if vz_sweep(j) > xsky6(4)
                z0_case(8) = max(0.75,zsw(8) - 0.25);
            else
                z0_case(8) = zsw(8) + 0.25;
            end

            % Scale switch guesses so all switches remain inside the new tf.
            scale_factor = z0_case(8)/zsw(8);
            z0_case(9:13) = zsw(9:13)*scale_factor;

            sens{j} = solve_switch_time_case(case_names{j},vz_sweep(j),xsky_case,x0,p,z0_case,opts_sw);
        end

        sj = sens{j};
        fprintf('\n%s\n', sj.caseName);
        fprintf('  target vz(tf)        = %.6f m/s\n', sj.vzTarget);
        fprintf('  exitflag             = %d\n', sj.exitflag);
        fprintf('  tf                   = %.6f s\n', sj.tf);
        fprintf('  T switches           = %.6f, %.6f s\n', sj.Tswitches(1), sj.Tswitches(2));
        fprintf('  tau switches         = %.6f, %.6f, %.6f s\n', sj.tauswitches(1), sj.tauswitches(2), sj.tauswitches(3));
        fprintf('  terminal residual    = %.6e\n', sj.terminalResidual);
        fprintf('  Hamiltonian final    = %.6e\n', sj.Hfinal);
        fprintf('  Hamiltonian variation= %.6e\n', sj.Hvariation);
    end

    % Summary table for the MATLAB command window.
    Case = strings(length(sens),1);
    vz_tf = zeros(length(sens),1);
    tf_s = zeros(length(sens),1);
    terminalResidual = zeros(length(sens),1);
    H_tf = zeros(length(sens),1);
    H_variation = zeros(length(sens),1);
    exitflag = zeros(length(sens),1);

    for j = 1:length(sens)
        sj = sens{j};
        Case(j) = string(sj.caseName);
        vz_tf(j) = sj.vzTarget;
        tf_s(j) = sj.tf;
        terminalResidual(j) = sj.terminalResidual;
        H_tf(j) = sj.Hfinal;
        H_variation(j) = sj.Hvariation;
        exitflag(j) = sj.exitflag;
    end

    SensitivityResults = table(Case,vz_tf,tf_s,terminalResidual,H_tf,H_variation,exitflag);
    disp(SensitivityResults);

    % Save table for easy use in the report.
    writetable(SensitivityResults,'phaseA_terminal_velocity_sensitivity.csv');

    % Trajectory comparison.
    figure;
    hold on;
    for j = 1:length(sens)
        sj = sens{j};
        plot(sj.x(:,1),sj.x(:,2),'LineWidth',1.8);
    end
    grid on; axis equal;
    xlabel('p_x [m]');
    ylabel('p_z [m]');
    title('Phase A Terminal Velocity Sensitivity: Trajectories');
    legend(case_names,'Location','best');

    % Physical thrust comparison.
    figure;
    subplot(2,1,1);
    hold on;
    for j = 1:length(sens)
        sj = sens{j};
        plot(sj.t,sj.Tphys,'LineWidth',1.5);
    end
    grid on;
    ylabel('T [N]');
    title('Terminal Velocity Sensitivity: Physical Thrust');
    legend(case_names,'Location','best');

    % Torque comparison.
    subplot(2,1,2);
    hold on;
    for j = 1:length(sens)
        sj = sens{j};
        plot(sj.t,sj.u(:,2),'LineWidth',1.5);
    end
    grid on;
    ylabel('\tau [N*m]');
    xlabel('Time [s]');
    title('Terminal Velocity Sensitivity: Torque');
    legend(case_names,'Location','best');

    % Hamiltonian comparison.
    figure;
    hold on;
    for j = 1:length(sens)
        sj = sens{j};
        plot(sj.t,sj.H,'LineWidth',1.5);
    end
    grid on;
    xlabel('Time [s]');
    ylabel('Hamiltonian');
    title('Terminal Velocity Sensitivity: Hamiltonian');
    legend(case_names,'Location','best');

end

%% ===== Local Functions =====

function A = makeA(g,m)
A = [0 0 1 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 g 0 0;
     0 0 0 0 0 0 -g/m;
     0 0 0 0 0 1 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0];
end

function B = makeB(m,I,alpha)
B = [0 0;
     0 0;
     0 0;
     1/m 0;
     0 0;
     0 1/I;
     -alpha 0];
end

function dy = canonical_ode(~,y,p)

x = y(1:7);
lam = y(8:14);

[u,~] = control_and_H(x,lam,p);

xdot = p.A*x + p.B*u;
lamdot = -p.A.'*lam;

dy = [xdot; lamdot];

end

function [u,H,phiT,phitau] = control_and_H(x,lam,p)

    phiT = (1/p.mh)*lam(4) - p.alpha*lam(7);
    phitau = (1/p.I)*lam(6);
    
    dT_center = 0.5*(p.dTmax+p.dTmin);
    dT_half   = 0.5*(p.dTmax-p.dTmin);
    
    tau_center = 0.5*(p.taumax+p.taumin);
    tau_half   = 0.5*(p.taumax-p.taumin);
    
    dT = dT_center - dT_half*tanh(phiT/p.eps);
    tau = tau_center - tau_half*tanh(phitau/p.eps);
    
    u = [dT; tau];
    
    H = 1 + lam.'*(p.A*x + p.B*u);
    
    %u = u.';
end

function r = phaseA_residual_massfree(z,x0,xsky6,p)

    lambda0 = z(1:7);
    tf = z(8);
    
    if tf <= 0
        r = 1e6*ones(8,1);
        return;
    end
    
    y0 = [x0; lambda0];
    
    [~,y] = ode45(@(t,y) canonical_ode(t,y,p), [0 tf], y0);
    
    yf = y(end,:).';
    xf = yf(1:7);
    lamf = yf(8:14);
    
    [u,Hf] = control_and_H(xf,lamf,p);
    
    r_unscaled = [
    xf(1:6) - xsky6;
    lamf(7);
    Hf
    ];

    scale = [
        20;    % px scale [m]
        25;    % pz scale [m]
        5;     % vx scale [m/s]
        5;     % vz scale [m/s]
        0.5;   % theta scale [rad]
        1.0;   % omega scale [rad/s]
        1.0;   % lambda_m scale
        1.0    % Hamiltonian scale
        ];
    
    r = r_unscaled ./ scale;

end

function r = switch_time_residual(z,x0,xsky6,p)

    lambda0 = z(1:7);
    tf      = z(8);
    tT      = z(9:10);
    ttau    = z(11:13);

    % Enforce switch-time ordering
    if tf <= 0 || ...
       tT(1) <= 0 || tT(2) <= tT(1) || tT(2) >= tf || ...
       ttau(1) <= 0 || ttau(2) <= ttau(1) || ...
       ttau(3) <= ttau(2) || ttau(3) >= tf
        r = 1e3*ones(13,1);
        return;
    end

    [t,y,u,H,phiT,phitau] = simulate_switch_solution(z,x0,p);

    xf   = y(end,1:7).';
    lamf = y(end,8:14).';

    % Switching functions at switch times
    phiT_sw1 = interp1(t,phiT,tT(1),'linear');
    phiT_sw2 = interp1(t,phiT,tT(2),'linear');

    phitau_sw1 = interp1(t,phitau,ttau(1),'linear');
    phitau_sw2 = interp1(t,phitau,ttau(2),'linear');
    phitau_sw3 = interp1(t,phitau,ttau(3),'linear');

    % Strong terminal accuracy
    terminal_scale = [0.1; 0.1; 0.25; 0.25; 0.025; 0.025];
    r_terminal = (xf(1:6)-xsky6)./terminal_scale;

    r_mass_free = lamf(7);
    r_H = H(end);

    % Switching residuals
    r_switch = [
        phiT_sw1/1e-1;
        phiT_sw2/1e-1;
        phitau_sw1/1e-1;
        phitau_sw2/1e-1;
        phitau_sw3/1e-1
    ];

    r = [
        r_terminal;
        r_mass_free;
        r_H;
        r_switch
    ];
end

function [t_all,y_all,u_all,H_all,phiT_all,phitau_all] = simulate_switch_solution(z,x0,p)

    lambda0 = z(1:7);
    tf      = z(8);
    tT      = z(9:10);
    ttau    = z(11:13);

    y0 = [x0; lambda0];

    % Include all switch times as integration breakpoints
    breaks = unique([0, tT(:).', ttau(:).', tf]);
    breaks = sort(breaks);

    t_all = [];
    y_all = [];

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

    yinit = y0;

    for k = 1:length(breaks)-1
        tspan = [breaks(k), breaks(k+1)];

        if abs(tspan(2)-tspan(1)) < 1e-10
            continue;
        end

        [tseg,yseg] = ode45(@(t,y) bangbang_canonical_ode(t,y,z,p), tspan, yinit, opts);

        if k > 1
            tseg = tseg(2:end);
            yseg = yseg(2:end,:);
        end

        t_all = [t_all; tseg];
        y_all = [y_all; yseg];

        yinit = yseg(end,:).';
    end

    u_all = zeros(length(t_all),2);
    H_all = zeros(length(t_all),1);
    phiT_all = zeros(length(t_all),1);
    phitau_all = zeros(length(t_all),1);

    for i = 1:length(t_all)
        x = y_all(i,1:7).';
        lam = y_all(i,8:14).';

        [ui,Hi,phiTi,phitaui] = bangbang_control_and_H(t_all(i),x,lam,z,p);

        u_all(i,:) = ui.';
        H_all(i) = Hi;
        phiT_all(i) = phiTi;
        phitau_all(i) = phitaui;
    end
end

function dy = bangbang_canonical_ode(t,y,z,p)

    x = y(1:7);
    lam = y(8:14);

    [u,~,~,~] = bangbang_control_and_H(t,x,lam,z,p);

    xdot = p.A*x + p.B*u;
    lamdot = -p.A.'*lam;

    dy = [xdot; lamdot];
end

function [u,H,phiT,phitau] = bangbang_control_and_H(t,x,lam,z,p)

    tf   = z(8);
    tT   = z(9:10);
    ttau = z(11:13);

    % Switching functions
    phiT = (1/p.mh)*lam(4) - p.alpha*lam(7);
    phitau = (1/p.I)*lam(6);

    % Thrust sequence: Tmax -> Tmin -> Tmax
    if t < tT(1)
        Tphys = p.Tmax;
    elseif t < tT(2)
        Tphys = p.Tmin;
    else
        Tphys = p.Tmax;
    end

    dT = Tphys - p.mh*p.g;

    % Torque sequence: taumax -> taumin -> taumax -> taumin
    if t < ttau(1)
        tau = p.taumax;
    elseif t < ttau(2)
        tau = p.taumin;
    elseif t < ttau(3)
        tau = p.taumax;
    else
        tau = p.taumin;
    end

    u = [dT; tau];

    H = 1 + lam.'*(p.A*x + p.B*u);
end

function result = solve_switch_time_case(caseName,vzTarget,xsky_case,x0,p,z0,opts_sw)

    lb = [-1e4*ones(7,1); 0.5; 0; 0; 0; 0; 0];
    ub = [ 1e4*ones(7,1); 20; 20; 20; 20; 20; 20];

    opts_case = optimoptions(opts_sw,'Display','off');

    resfun = @(z) switch_time_residual(z,x0,xsky_case,p);

    [z,resnorm,resid,exitflag,output] = lsqnonlin(resfun,z0,lb,ub,opts_case);

    result = package_switch_result(caseName,vzTarget,z,exitflag,resid,xsky_case,x0,p);
    result.resnorm_lsq = resnorm;
    result.output = output;

    if exitflag <= 0 || result.terminalResidual > 1e-6 || abs(result.Hfinal) > 1e-5
        warning(['Sensitivity case "',caseName,'" may not be a valid exact solution. ', ...
                 'Check residuals and consider changing the switching structure.']);
    end
end

function result = package_switch_result(caseName,vzTarget,z,exitflag,resid,xsky_case,x0,p)

    [t,y,u,H,phiT,phitau] = simulate_switch_solution(z,x0,p);

    x = y(:,1:7);
    lam = y(:,8:14);
    Tphys = p.mh*p.g + u(:,1);

    terminalError = x(end,1:6).' - xsky_case;

    result.caseName = caseName;
    result.vzTarget = vzTarget;
    result.z = z;
    result.exitflag = exitflag;
    result.scaledResidualNorm = norm(resid);
    result.terminalResidual = norm(terminalError);
    result.terminalError = terminalError;
    result.tf = z(8);
    result.Tswitches = z(9:10);
    result.tauswitches = z(11:13);
    result.Hfinal = H(end);
    result.Hinitial = H(1);
    result.Hvariation = max(H)-min(H);
    result.maxT = max(Tphys);
    result.minT = min(Tphys);
    result.maxAbsTau = max(abs(u(:,2)));
    result.t = t;
    result.y = y;
    result.x = x;
    result.lam = lam;
    result.u = u;
    result.Tphys = Tphys;
    result.H = H;
    result.phiT = phiT;
    result.phitau = phitau;

    % Default fields used by the sensitivity-study structure array.
    % These keep every entry in sens(j) structurally identical, including
    % the baseline case that reuses the already-computed solution.
    result.resnorm_lsq = NaN;
    result.output = struct();
end
