%% Part IV: Nonlinear Optimal Control Near Hover
% Compares:
% 1) Linear optimal controller on linear dynamics
% 2) Same linear controller on nonlinear dynamics
% 3) Nonlinear optimal controller on nonlinear dynamics

clear; clc; close all;

%% =========================
% Parameters
% =========================

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

p.tf = 2.5;
p.n = 7;

p.xhover = [0; 50; 0; 0; 0; 0; p.mh];

p.Ah = makeA_hover(p.g,p.mh);
p.Bh = makeB_hover(p.mh,p.I,p.alpha);

% Use tuned feasible weights
p.Q = diag([10, 10, 1, 1, 50, 5, 0.1]);
p.R = diag([0.05, 2.7]);

% No terminal cost in Part IV
p.Qf = zeros(7);

% Near-hover initial condition
e0 = [1.5; -1; 0.3; -0.2; deg2rad(2); deg2rad(0.5); 0];
x0 = p.xhover + e0;

tgrid = linspace(0,p.tf,1001);

%% =========================
% Case 1: Linear LQR on Linear Dynamics
% =========================

linRes = solve_linear_lqr(p,e0,tgrid);

%% =========================
% Case 2: Same Linear Controller on Nonlinear Dynamics
% =========================

nlLinRes = simulate_nonlinear_with_linear_controller(p,x0,linRes,tgrid);

%% =========================
% Case 3: Nonlinear Optimal Control by Indirect Shooting
% =========================

nlOptRes = solve_nonlinear_shooting(p,x0,e0,linRes,tgrid);

%% =========================
% Print Comparison
% =========================

print_partIV_summary(linRes,nlLinRes,nlOptRes,p);

%% =========================
% Plot Comparison
% =========================

plot_partIV_results(linRes,nlLinRes,nlOptRes,p);

%% ============================================================
% Local Functions
% ============================================================

function A = makeA_hover(g,m)
A = [0 0 1 0 0 0 0;
     0 0 0 1 0 0 0;
     0 0 0 0 g 0 0;
     0 0 0 0 0 0 -g/m;
     0 0 0 0 0 1 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0];
end

function B = makeB_hover(m,I,alpha)
B = [0 0;
     0 0;
     0 0;
     1/m 0;
     0 0;
     0 1/I;
     -alpha 0];
end

function res = solve_linear_lqr(p,e0,tgrid)

    Ptf = p.Qf;
    pvec0 = Ptf(:);

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

    [tPraw,Praw] = ode45(@(t,Pvec) riccati_ode(t,Pvec,p), [p.tf 0], pvec0, opts);

    [tP,idx] = sort(tPraw);
    Pvec = Praw(idx,:);

    [tx,eraw] = ode45(@(t,e) linear_closed_loop_ode(t,e,p,tP,Pvec), [0 p.tf], e0, opts);

    e = interp1(tx,eraw,tgrid,'pchip').';

    u = zeros(2,length(tgrid));
    H = zeros(1,length(tgrid));

    for k = 1:length(tgrid)
        Pk = getP(tgrid(k),tP,Pvec,p.n);
        K = p.R\(p.Bh.'*Pk);
        u(:,k) = -K*e(:,k);

        lambda = 2*Pk*e(:,k);
        edot = p.Ah*e(:,k)+p.Bh*u(:,k);
        H(k) = e(:,k).'*p.Q*e(:,k) + u(:,k).'*p.R*u(:,k) + lambda.'*edot;
    end

    res.name = 'Linear on Linear';
    res.t = tgrid;
    res.e = e;
    res.x = p.xhover + e;
    res.u = u;
    res.H = H;
    res.P_t = tP;
    res.P_vec = Pvec;
end

function dPvec = riccati_ode(t,Pvec,p)

    P = reshape(Pvec,p.n,p.n);
    P = 0.5*(P+P.');

    A = p.Ah;
    B = p.Bh;

    Pdot = -(A.'*P + P*A - P*B*(p.R\B.')*P + p.Q);

    dPvec = Pdot(:);
end

function de = linear_closed_loop_ode(t,e,p,tP,Pvec)

    P = getP(t,tP,Pvec,p.n);
    K = p.R\(p.Bh.'*P);
    u = -K*e;

    de = p.Ah*e + p.Bh*u;
end

function res = simulate_nonlinear_with_linear_controller(p,x0,linRes,tgrid)

    opts = odeset('RelTol',1e-9,'AbsTol',1e-11);

    rhs = @(t,x) nonlinear_with_linear_control_ode(t,x,p,linRes);

    [tx,xraw] = ode45(rhs,[0 p.tf],x0,opts);

    x = interp1(tx,xraw,tgrid,'pchip').';
    e = x - p.xhover;

    u = zeros(2,length(tgrid));
    H = zeros(1,length(tgrid));

    for k = 1:length(tgrid)
        t = tgrid(k);
        ek = e(:,k);
        Pk = getP(t,linRes.P_t,linRes.P_vec,p.n);
        K = p.R\(p.Bh.'*Pk);
        u(:,k) = -K*ek;

        lambda = 2*Pk*ek;
        f = nonlinear_f(x(:,k),u(:,k),p);
        H(k) = ek.'*p.Q*ek + u(:,k).'*p.R*u(:,k) + lambda.'*f;
    end

    res.name = 'Linear Controller on Nonlinear';
    res.t = tgrid;
    res.x = x;
    res.e = e;
    res.u = u;
    res.H = H;
end

function dx = nonlinear_with_linear_control_ode(t,x,p,linRes)

    e = x - p.xhover;

    P = getP(t,linRes.P_t,linRes.P_vec,p.n);
    K = p.R\(p.Bh.'*P);
    u = -K*e;

    u = saturate_u(u,p);

    dx = nonlinear_f(x,u,p);
end

function res = solve_nonlinear_shooting(p,x0,e0,linRes,tgrid)

    % Initial costate guess from linear Riccati solution
    P0 = getP(0,linRes.P_t,linRes.P_vec,p.n);
    lambda0_guess = 2*P0*e0;

    lb = -1e4*ones(7,1);
    ub =  1e4*ones(7,1);

    opts = optimoptions('lsqnonlin', ...
        'Display','iter', ...
        'MaxIterations',1000, ...
        'MaxFunctionEvaluations',2e4, ...
        'FunctionTolerance',1e-10, ...
        'StepTolerance',1e-10, ...
        'OptimalityTolerance',1e-10);

    residual_fun = @(lambda0) nonlinear_shooting_residual(lambda0,x0,p);

    [lambda0_sol,resnorm,resid,exitflag,output] = lsqnonlin(residual_fun,lambda0_guess,lb,ub,opts);

    fprintf('\nNonlinear shooting exitflag = %d\n',exitflag);
    fprintf('Nonlinear terminal costate residual norm = %.6e\n',norm(resid));

    optsODE = odeset('RelTol',1e-9,'AbsTol',1e-11);

    y0 = [x0; lambda0_sol];

    [tx,yraw] = ode45(@(t,y) nonlinear_canonical_ode(t,y,p), [0 p.tf], y0, optsODE);

    y = interp1(tx,yraw,tgrid,'pchip').';

    x = y(1:7,:);
    lambda = y(8:14,:);
    e = x - p.xhover;

    u = zeros(2,length(tgrid));
    H = zeros(1,length(tgrid));
    statRes = zeros(2,length(tgrid));

    for k = 1:length(tgrid)
        [u(:,k),H(k),statRes(:,k)] = nonlinear_control_H_stationarity(x(:,k),lambda(:,k),p);
    end

    res.name = 'Nonlinear Optimal on Nonlinear';
    res.t = tgrid;
    res.x = x;
    res.e = e;
    res.lambda = lambda;
    res.u = u;
    res.H = H;
    res.stationarity = statRes;
    res.terminalResidual = lambda(:,end);
end

function r = nonlinear_shooting_residual(lambda0,x0,p)

    optsODE = odeset('RelTol',1e-8,'AbsTol',1e-10);

    y0 = [x0; lambda0];

    [~,y] = ode45(@(t,y) nonlinear_canonical_ode(t,y,p), [0 p.tf], y0, optsODE);

    yf = y(end,:).';

    lambda_tf = yf(8:14);

    % Terminal condition: lambda(tf)=0
    r = lambda_tf;
end

function dy = nonlinear_canonical_ode(~,y,p)

    x = y(1:7);
    lambda = y(8:14);

    [u,~,~] = nonlinear_control_H_stationarity(x,lambda,p);

    dx = nonlinear_f(x,u,p);

    e = x - p.xhover;

    T = p.mh*p.g + u(1);
    m = x(7);
    th = x(5);

    lambdadot = zeros(7,1);

    lambdadot(1) = -2*p.Q(1,1)*e(1);
    lambdadot(2) = -2*p.Q(2,2)*e(2);
    lambdadot(3) = -(2*p.Q(3,3)*e(3) + lambda(1));
    lambdadot(4) = -(2*p.Q(4,4)*e(4) + lambda(2));

    lambdadot(5) = -2*p.Q(5,5)*e(5) ...
                   - lambda(3)*(T/m)*cos(th) ...
                   + lambda(4)*(T/m)*sin(th);

    lambdadot(6) = -(2*p.Q(6,6)*e(6) + lambda(5));

    lambdadot(7) = -2*p.Q(7,7)*e(7) ...
                   + lambda(3)*T*sin(th)/(m^2) ...
                   + lambda(4)*T*cos(th)/(m^2);

    dy = [dx; lambdadot];
end

function dx = nonlinear_f(x,u,p)

    px = x(1);
    pz = x(2);
    vx = x(3);
    vz = x(4);
    th = x(5);
    om = x(6);
    m  = x(7);

    dT  = u(1);
    tau = u(2);

    T = p.mh*p.g + dT;

    dx = zeros(7,1);

    dx(1) = vx;
    dx(2) = vz;
    dx(3) = (T/m)*sin(th);
    dx(4) = (T/m)*cos(th) - p.g;
    dx(5) = om;
    dx(6) = tau/p.I;
    dx(7) = -p.alpha*T;
end

function [u,H,statRes] = nonlinear_control_H_stationarity(x,lambda,p)

    th = x(5);
    m = x(7);

    RT = p.R(1,1);
    Rtau = p.R(2,2);

    dT = -(lambda(3)*sin(th)/m + lambda(4)*cos(th)/m - p.alpha*lambda(7))/(2*RT);
    tau = -lambda(6)/(2*p.I*Rtau);

    u = [dT; tau];

    u = saturate_u(u,p);

    e = x - p.xhover;
    f = nonlinear_f(x,u,p);

    H = e.'*p.Q*e + u.'*p.R*u + lambda.'*f;

    statRes = zeros(2,1);
    statRes(1) = 2*RT*u(1) + lambda(3)*sin(th)/m + lambda(4)*cos(th)/m - p.alpha*lambda(7);
    statRes(2) = 2*Rtau*u(2) + lambda(6)/p.I;
end

function u = saturate_u(u,p)

    u(1) = min(max(u(1),p.dTmin),p.dTmax);
    u(2) = min(max(u(2),p.taumin),p.taumax);
end

function P = getP(t,tP,Pvec,n)

    pinterp = interp1(tP,Pvec,t,'pchip');
    P = reshape(pinterp.',n,n);
    P = 0.5*(P+P.');
end

function print_partIV_summary(linRes,nlLinRes,nlOptRes,p)

    results = {linRes,nlLinRes,nlOptRes};

    fprintf('\n================ PART IV COMPARISON ================\n');
    fprintf('%-32s %-12s %-12s %-12s %-12s %-12s\n', ...
        'Case','J','||e(tf)||','maxT','maxTau','Hvar');

    for i = 1:length(results)
        res = results{i};

        Tphys = p.mh*p.g + res.u(1,:);

        integrand = zeros(1,length(res.t));
        for k = 1:length(res.t)
            ek = res.e(:,k);
            uk = res.u(:,k);
            integrand(k) = ek.'*p.Q*ek + uk.'*p.R*uk;
        end

        J = trapz(res.t,integrand);

        fprintf('%-32s %-12.4e %-12.4e %-12.4e %-12.4e %-12.4e\n', ...
            res.name,J,norm(res.e(:,end)),max(Tphys),max(abs(res.u(2,:))),max(res.H)-min(res.H));
    end

    fprintf('\nNonlinear optimal terminal costate residual: %.6e\n',norm(nlOptRes.terminalResidual));
    fprintf('Nonlinear optimal stationarity residual: %.6e\n',max(vecnorm(nlOptRes.stationarity,2,1)));
end

function plot_partIV_results(linRes,nlLinRes,nlOptRes,p)

    results = {linRes,nlLinRes,nlOptRes};
    names = {linRes.name,nlLinRes.name,nlOptRes.name};

    labels = {'p_x error','p_z error','v_x','v_z','theta','omega','mass error'};

    figure('Name','Part IV State Error Comparison');
    tiledlayout(4,2);

    for j = 1:7
        nexttile; hold on; grid on;
        for i = 1:length(results)
            plot(results{i}.t,results{i}.e(j,:),'LineWidth',1.4);
        end
        xlabel('Time [s]');
        ylabel(labels{j});
        title(labels{j});
        legend(names,'Location','best');
    end

    figure('Name','Part IV Control Comparison');
    tiledlayout(2,1);

    nexttile; hold on; grid on;
    for i = 1:length(results)
        Tphys = p.mh*p.g + results{i}.u(1,:);
        plot(results{i}.t,Tphys,'LineWidth',1.5);
    end
    yline(p.Tmin,'--');
    yline(p.Tmax,'--');
    xlabel('Time [s]');
    ylabel('T [N]');
    title('Physical Thrust Comparison');
    legend(names,'Location','best');

    nexttile; hold on; grid on;
    for i = 1:length(results)
        plot(results{i}.t,results{i}.u(2,:),'LineWidth',1.5);
    end
    yline(p.taumin,'--');
    yline(p.taumax,'--');
    xlabel('Time [s]');
    ylabel('\tau [N*m]');
    title('Torque Comparison');
    legend(names,'Location','best');

    figure('Name','Part IV Trajectory Comparison');
    hold on; grid on; axis equal;
    for i = 1:length(results)
        plot(results{i}.x(1,:),results{i}.x(2,:),'LineWidth',1.8);
    end
    xlabel('p_x [m]');
    ylabel('p_z [m]');
    title('Position Trajectory Comparison');
    legend(names,'Location','best');

    figure('Name','Part IV Hamiltonian Comparison');
    hold on; grid on;
    for i = 1:length(results)
        plot(results{i}.t,results{i}.H,'LineWidth',1.5);
    end
    xlabel('Time [s]');
    ylabel('H');
    title('Hamiltonian Evolution');
    legend(names,'Location','best');
end