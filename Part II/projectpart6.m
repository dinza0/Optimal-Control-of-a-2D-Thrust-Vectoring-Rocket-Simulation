%% Part II: Optimal Trajectory Tracking
% Uses a Part I optimal state trajectory as x_ref(t)
% Then designs an affine finite-horizon LQR tracking controller.

clear; clc; close all;

%% =========================
%  1. Configuration
%  =========================

cfg.g      = 9.81;
cfg.I      = 10;
cfg.m0     = 20;
cfg.alpha  = 0.005;

cfg.Tmin   = 0;
cfg.Tmax   = 300;
cfg.taumin = -30;
cfg.taumax = 30;
cfg.mdry   = 5;

cfg.tf     = 5.0;
cfg.n      = 7;
cfg.r      = 2;
cfg.Nplot  = 1001;

tgrid = linspace(0, cfg.tf, cfg.Nplot);

% Weak descent nominal trajectory
cfg.weak_px0 = 0;
cfg.weak_pz0 = 50;
cfg.weak_vx  = 10;
cfg.weak_vz  = -3;

%% =========================
%  2. Tuned Part I Weights
%  =========================

Q  = diag([10, 10, 1, 1, 50, 5, 0.1]);
R  = diag([0.04, 7.5]);
Qf = diag([100, 100, 10, 10, 100, 20, 1]);

%% =========================
%  3. Define Weak LTV Model
%  =========================

model.name = 'Weak LTV Descent Tracking';
model.Afun = @(t) makeA(cfg.g, mstar_weak(t, cfg));
model.Bfun = @(t) makeB(mstar_weak(t, cfg), cfg.I, cfg.alpha);
model.refPhysical = @(t) reference_weak_physical(t, cfg);

%% =========================
%  4. Generate Part I Reference Trajectory
%  =========================
% This is the state trajectory generated in Part I.
% It is a perturbation trajectory, not the absolute physical trajectory.

x0_partI = [8; -4; 2; -1; deg2rad(5); deg2rad(2); 0];

partI_ref = solve_partI_lqr_reference(model, cfg, Q, R, Qf, x0_partI, tgrid);

xref_grid = partI_ref.x;

fprintf('\nGenerated Part I reference trajectory.\n');
fprintf('||xref(0)||  = %.6e\n', norm(xref_grid(:,1)));
fprintf('||xref(tf)|| = %.6e\n', norm(xref_grid(:,end)));

%% =========================
%  5. Choose Different Tracking Initial Condition
%  =========================
% This should be substantially different from xref(0).

x0_track = [-6; 6; -3; 2; deg2rad(-8); deg2rad(-3); 0];

fprintf('\nTracking initial condition chosen.\n');
fprintf('||x0_track - xref(0)|| = %.6e\n', norm(x0_track - xref_grid(:,1)));

%% =========================
%  6. Solve Tracking Controller
%  =========================

trackRes = solve_tracking_lqr(model, cfg, Q, R, xref_grid, x0_track, tgrid);

print_tracking_summary(trackRes, cfg);

%% =========================
%  7. Run Sensitivity Cases
%  =========================

icCases = {};
icCases{end+1}.name = 'Nominal tracking IC';
icCases{end}.x0 = x0_track;

icCases{end+1}.name = '25 percent larger IC';
icCases{end}.x0 = 1.25*x0_track;

icCases{end+1}.name = '50 percent larger IC';
icCases{end}.x0 = 1.50*x0_track;

icCases{end+1}.name = 'High vx tracking error';
icCases{end}.x0 = [-6; 6; -6; 2; deg2rad(-8); deg2rad(-3); 0];

icCases{end+1}.name = 'High theta tracking error';
icCases{end}.x0 = [-6; 6; -3; 2; deg2rad(-15); deg2rad(-3); 0];

summaryRows = {};

for k = 1:numel(icCases)
    res = solve_tracking_lqr(model, cfg, Q, R, xref_grid, icCases{k}.x0, tgrid);

    feasible = ~(res.thrust_violation || res.torque_violation || ...
                 res.altitude_violation || res.mass_violation);

    summaryRows = [summaryRows;
        {icCases{k}.name, norm(icCases{k}.x0 - xref_grid(:,1)), ...
         norm(res.err(:,end)), max(vecnorm(res.err,2,1)), ...
         res.Jtotal, res.Jcontrol, max(res.Tphys), ...
         max(abs(res.tauphys)), feasible}];
end

sensitivityTable = cell2table(summaryRows, ...
    'VariableNames', {'IC_Case','InitialTrackingError','FinalTrackingError', ...
                      'PeakTrackingError','J_total','J_control', ...
                      'MaxThrust','MaxAbsTorque','Feasible'});

disp(sensitivityTable);
writetable(sensitivityTable, 'partII_tracking_sensitivity_table.csv');

%% =========================
%  8. Plot Results
%  =========================

plot_tracking_results(trackRes, xref_grid, cfg, model);

%% ============================================================
%  Local Functions
%  ============================================================

function res = solve_partI_lqr_reference(model, cfg, Q, R, Qf, x0, tgrid)

    n = cfg.n;

    Ptf = Qf;
    p0_vec = Ptf(:);

    rhs = @(t,pvec) riccati_regulation_ode(t, pvec, model.Afun, model.Bfun, Q, R, n);

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    [tP_raw, Pvec_raw] = ode45(rhs, [cfg.tf 0], p0_vec, opts);

    [tP, idx] = sort(tP_raw);
    Pvec = Pvec_raw(idx,:);

    cl_rhs = @(t,x) closed_loop_regulation(t, x, model.Afun, model.Bfun, R, tP, Pvec, n);

    [tx, xraw] = ode45(cl_rhs, [0 cfg.tf], x0, opts);
    xgrid = interp1(tx, xraw, tgrid, 'pchip').';

    ugrid = zeros(cfg.r, numel(tgrid));
    for k = 1:numel(tgrid)
        t = tgrid(k);
        B = model.Bfun(t);
        P = getP(t, tP, Pvec, n);
        K = R \ (B.'*P);
        ugrid(:,k) = -K*xgrid(:,k);
    end

    res.t = tgrid;
    res.x = xgrid;
    res.u = ugrid;
end

function res = solve_tracking_lqr(model, cfg, Q, R, xref_grid, x0, tgrid)

    n = cfg.n;

    % No terminal cost in the Part II tracking problem statement:
    Ptf = zeros(n,n);
    stf = zeros(n,1);

    ztf = [Ptf(:); stf];

    rhs = @(t,z) tracking_riccati_ode(t, z, model.Afun, model.Bfun, Q, R, xref_grid, tgrid, n);

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
    [tz_raw, zraw] = ode45(rhs, [cfg.tf 0], ztf, opts);

    [tz, idx] = sort(tz_raw);
    zsol = zraw(idx,:);

    cl_rhs = @(t,x) tracking_closed_loop_ode(t, x, model.Afun, model.Bfun, R, tz, zsol, n);

    [tx, xraw] = ode45(cl_rhs, [0 cfg.tf], x0, opts);
    xgrid = interp1(tx, xraw, tgrid, 'pchip').';

    ugrid = zeros(cfg.r, numel(tgrid));
    lambdagrid = zeros(n, numel(tgrid));
    errgrid = zeros(n, numel(tgrid));
    Hgrid = zeros(1, numel(tgrid));
    stationarity = zeros(cfg.r, numel(tgrid));

    xabs_grid = zeros(n, numel(tgrid));
    Tphys = zeros(1, numel(tgrid));
    tauphys = zeros(1, numel(tgrid));

    state_integrand = zeros(1, numel(tgrid));
    control_integrand = zeros(1, numel(tgrid));

    for k = 1:numel(tgrid)
        t = tgrid(k);
        x = xgrid(:,k);
        xref = xref_grid(:,k);

        [P, s] = getTrackingPS(t, tz, zsol, n);

        A = model.Afun(t);
        B = model.Bfun(t);

        u = -R \ (B.'*(P*x + s));
        lambda = 2*(P*x + s);

        err = x - xref;
        xdot = A*x + B*u;

        H = err.'*Q*err + u.'*R*u + lambda.'*xdot;

        ru = 2*R*u + B.'*lambda;

        physicalRef = model.refPhysical(t);
        xabs = physicalRef.x + x;

        Tnom = physicalRef.Tnom;
        Tphys(k) = Tnom + u(1);
        tauphys(k) = u(2);

        ugrid(:,k) = u;
        lambdagrid(:,k) = lambda;
        errgrid(:,k) = err;
        Hgrid(k) = H;
        stationarity(:,k) = ru;
        xabs_grid(:,k) = xabs;

        state_integrand(k) = err.'*Q*err;
        control_integrand(k) = u.'*R*u;
    end

    Jstate = trapz(tgrid, state_integrand);
    Jcontrol = trapz(tgrid, control_integrand);
    Jtotal = Jstate + Jcontrol;

    terminal_costate_res = lambdagrid(:,end);

    res.t = tgrid;
    res.x = xgrid;
    res.xref = xref_grid;
    res.err = errgrid;
    res.u = ugrid;
    res.lambda = lambdagrid;
    res.H = Hgrid;
    res.xabs = xabs_grid;
    res.Tphys = Tphys;
    res.tauphys = tauphys;

    res.Jstate = Jstate;
    res.Jcontrol = Jcontrol;
    res.Jtotal = Jtotal;

    res.maxStationarity = max(vecnorm(stationarity,2,1));
    res.terminalCostateNorm = norm(terminal_costate_res);

    res.maxTrackingError = max(vecnorm(errgrid,2,1));
    res.finalTrackingError = norm(errgrid(:,end));
    res.energy_uTu = trapz(tgrid, sum(ugrid.^2,1));

    res.minTphys = min(Tphys);
    res.maxTphys = max(Tphys);
    res.minTau = min(tauphys);
    res.maxTau = max(tauphys);

    res.minPz = min(xabs_grid(2,:));
    res.minMass = min(xabs_grid(7,:));

    res.thrust_violation = min(Tphys) < cfg.Tmin || max(Tphys) > cfg.Tmax;
    res.torque_violation = min(tauphys) < cfg.taumin || max(tauphys) > cfg.taumax;
    res.altitude_violation = res.minPz < 0;
    res.mass_violation = res.minMass < cfg.mdry;
end

function dz = tracking_riccati_ode(t, z, Afun, Bfun, Q, R, xref_grid, tgrid, n)

    P = reshape(z(1:n*n), n, n);
    P = 0.5*(P+P.');
    s = z(n*n+1:end);

    A = Afun(t);
    B = Bfun(t);
    G = B*(R\B.');

    xref = interp1(tgrid, xref_grid.', t, 'pchip').';

    % Riccati equation:
    % -Pdot = A'P + PA - P G P + Q
    Pdot = -(A.'*P + P*A - P*G*P + Q);

    % Tracking vector equation:
    % sdot = -A's + P G s + Q xref
    sdot = -A.'*s + P*G*s + Q*xref;

    dz = [Pdot(:); sdot];
end

function dx = tracking_closed_loop_ode(t, x, Afun, Bfun, R, tz, zsol, n)

    [P, s] = getTrackingPS(t, tz, zsol, n);

    A = Afun(t);
    B = Bfun(t);

    u = -R \ (B.'*(P*x + s));

    dx = A*x + B*u;
end

function dp = riccati_regulation_ode(t, pvec, Afun, Bfun, Q, R, n)

    P = reshape(pvec, n, n);
    P = 0.5*(P+P.');

    A = Afun(t);
    B = Bfun(t);

    Pdot = -(A.'*P + P*A - P*B*(R\B.')*P + Q);

    dp = Pdot(:);
end

function dx = closed_loop_regulation(t, x, Afun, Bfun, R, tP, Pvec, n)

    A = Afun(t);
    B = Bfun(t);
    P = getP(t, tP, Pvec, n);

    u = -R \ (B.'*P*x);

    dx = A*x + B*u;
end

function P = getP(t, tP, Pvec, n)

    p = interp1(tP, Pvec, t, 'pchip');
    P = reshape(p.', n, n);
    P = 0.5*(P+P.');
end

function [P, s] = getTrackingPS(t, tz, zsol, n)

    z = interp1(tz, zsol, t, 'pchip').';

    P = reshape(z(1:n*n), n, n);
    P = 0.5*(P+P.');
    s = z(n*n+1:end);
end

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

function m = mstar_weak(t, cfg)

    m = cfg.m0 * exp(-cfg.alpha*cfg.g*t);
end

function ref = reference_weak_physical(t, cfg)

    m = mstar_weak(t, cfg);

    ref.x = [cfg.weak_px0 + cfg.weak_vx*t;
             cfg.weak_pz0 + cfg.weak_vz*t;
             cfg.weak_vx;
             cfg.weak_vz;
             0;
             0;
             m];

    ref.Tnom = m*cfg.g;
end

function print_tracking_summary(res, cfg)

    fprintf('\n================ TRACKING SUMMARY ================\n');
    fprintf('J_state                      = %.6e\n', res.Jstate);
    fprintf('J_control                    = %.6e\n', res.Jcontrol);
    fprintf('J_total                      = %.6e\n', res.Jtotal);
    fprintf('Max tracking error norm       = %.6e\n', res.maxTrackingError);
    fprintf('Final tracking error norm     = %.6e\n', res.finalTrackingError);
    fprintf('Integral u^T u dt             = %.6e\n', res.energy_uTu);

    fprintf('\n--- Control Bounds ---\n');
    fprintf('min physical T                = %.6e N\n', res.minTphys);
    fprintf('max physical T                = %.6e N\n', res.maxTphys);
    fprintf('min tau                       = %.6e N*m\n', res.minTau);
    fprintf('max tau                       = %.6e N*m\n', res.maxTau);

    if res.thrust_violation
        fprintf('WARNING: thrust bound violated. Bounds: [%.2f, %.2f] N\n', cfg.Tmin, cfg.Tmax);
    else
        fprintf('Thrust bounds satisfied.\n');
    end

    if res.torque_violation
        fprintf('WARNING: torque bound violated. Bounds: [%.2f, %.2f] N*m\n', cfg.taumin, cfg.taumax);
    else
        fprintf('Torque bounds satisfied.\n');
    end

    fprintf('\n--- State Constraints ---\n');
    fprintf('min altitude pz               = %.6e m\n', res.minPz);
    fprintf('min mass m                    = %.6e kg\n', res.minMass);

    if res.altitude_violation
        fprintf('WARNING: altitude violated.\n');
    else
        fprintf('Altitude constraint satisfied.\n');
    end

    if res.mass_violation
        fprintf('WARNING: dry mass violated.\n');
    else
        fprintf('Dry mass constraint satisfied.\n');
    end

    fprintf('\n--- PMP Validation ---\n');
    fprintf('max ||2Ru + B^T lambda||      = %.6e\n', res.maxStationarity);
    fprintf('||lambda(tf)||                = %.6e\n', res.terminalCostateNorm);
    fprintf('Hamiltonian initial           = %.6e\n', res.H(1));
    fprintf('Hamiltonian final             = %.6e\n', res.H(end));
    fprintf('Hamiltonian variation         = %.6e\n', max(res.H)-min(res.H));
end

function plot_tracking_results(res, xref_grid, cfg, model)

    t = res.t;

    labels = {'p_x','p_z','v_x','v_z','\theta','\omega','m'};

    figure('Name','Part II Tracking Error States');
    tiledlayout(4,2);

    for k = 1:7
        nexttile;
        plot(t, res.err(k,:), 'LineWidth', 1.5);
        grid on;
        xlabel('Time [s]');
        ylabel(['e_' labels{k}]);
        title(['Tracking Error: e_' labels{k}]);
    end

    figure('Name','Part II State vs Reference');
    tiledlayout(4,2);

    for k = 1:7
        nexttile;
        plot(t, res.x(k,:), 'LineWidth', 1.5); hold on;
        plot(t, xref_grid(k,:), '--', 'LineWidth', 1.5);
        grid on;
        xlabel('Time [s]');
        ylabel(labels{k});
        title([labels{k} ': actual vs reference']);
        legend('Tracked state','Reference');
    end

    figure('Name','Part II Controls');
    tiledlayout(3,1);

    nexttile;
    plot(t, res.u(1,:), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('T perturbation [N]');
    title('Tracking Thrust Perturbation');

    nexttile;
    plot(t, res.u(2,:), 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('\tau [N*m]');
    title('Tracking Torque');

    nexttile;
    plot(t, res.Tphys, 'LineWidth', 1.5); hold on;
    yline(cfg.Tmin,'--');
    yline(cfg.Tmax,'--');
    grid on;
    xlabel('Time [s]');
    ylabel('T_{phys} [N]');
    title('Physical Thrust with Bounds');

    figure('Name','Part II Tracking Error Norm');
    plot(t, vecnorm(res.err,2,1), 'LineWidth', 1.7);
    grid on;
    xlabel('Time [s]');
    ylabel('||x-x_{ref}||_2');
    title('Tracking Error Norm');

    figure('Name','Part II Hamiltonian');
    plot(t, res.H, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('Hamiltonian');
    title('Tracking Hamiltonian Evolution');

    figure('Name','Part II Physical Mission Trajectory');
    plot(res.xabs(1,:), res.xabs(2,:), 'LineWidth', 1.8); hold on;
    grid on;
    xlabel('p_x [m]');
    ylabel('p_z [m]');
    title('Tracking Physical Mission Trajectory with Attitude Vectors');
    axis equal;

    stride = max(1, floor(numel(t)/20));
    scale = 2.0;

    for k = 1:stride:numel(t)
        px = res.xabs(1,k);
        pz = res.xabs(2,k);
        th = res.xabs(5,k);

        dx = scale*sin(th);
        dz = scale*cos(th);

        quiver(px, pz, dx, dz, 0, 'LineWidth', 1.0, 'MaxHeadSize', 1.5);
    end
end