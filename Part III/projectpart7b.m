%% Part III Phase B: Return-and-Land Energy Optimal Problem
% Solves both landing manifolds using PMP + bvp4c

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

% Energy weights
p.R = diag([0.05, 2.7]);
p.rho = 1.0;

% Circular platform parameters
p.pc = 0;
p.r  = 2;

%% Initial state for Phase B
% Use your Phase A sky rendezvous result.
% If Phase A target was xsky6 = [10;15;0;0;0;0], use:
dm_sky = 0; % replace with xbb(end,7) from Phase A if available

x_sky = [10; 25; 0; 0; 0; 0; dm_sky];

%% Solve landing manifold 1
fprintf('\n=============================\n');
fprintf('Solving Phase B: Manifold M1\n');
fprintf('=============================\n');

p.land_px = 0;   % desired ground landing x-position
resM1 = solve_phaseB_bvp('M1x',x_sky,p);

%% Solve landing manifold 2
fprintf('\n=============================\n');
fprintf('Solving Phase B: Manifold M2\n');
fprintf('=============================\n');

resM2 = solve_phaseB_bvp('M2',x_sky,p);

%% Compare
fprintf('\n\n============= PHASE B COMPARISON =============\n');
fprintf('%-10s %-12s %-12s %-12s %-12s %-12s\n', ...
    'Manifold','tf','J','maxT','maxTau','H(tf)');

fprintf('%-10s %-12.5f %-12.5e %-12.5f %-12.5f %-12.3e\n', ...
    'M1',resM1.tf,resM1.J,max(resM1.Tphys),max(abs(resM1.u(:,2))),resM1.H(end));

fprintf('%-10s %-12.5f %-12.5e %-12.5f %-12.5f %-12.3e\n', ...
    'M2',resM2.tf,resM2.J,max(resM2.Tphys),max(abs(resM2.u(:,2))),resM2.H(end));

%% Plot comparison
figure;
plot(resM1.x(:,1),resM1.x(:,2),'LineWidth',2); hold on;
plot(resM2.x(:,1),resM2.x(:,2),'LineWidth',2);

% draw circular landing platform
ang = linspace(0,2*pi,200);
xc = p.pc + p.r*cos(ang);
zc = 0    + p.r*sin(ang);
plot(xc,zc,'r--','LineWidth',2);
patch(xc,zc,[1 0.85 0.85],'FaceAlpha',0.20,'EdgeColor','r');

grid on; axis equal;
xlabel('p_x [m]');
ylabel('p_z [m]');
legend('M1 Ground landing','M2 Circular platform trajectory','Landing platform');
title('Phase B Landing Trajectories');

%% ============================================================
% Local functions
% ============================================================

function res = solve_phaseB_bvp(manifold,x_sky,p)

    % Initial guess
    tf_guess = 12;

    s_mesh = linspace(0,1,60);

    switch manifold
        case 'M1x'
            xf_guess = [p.land_px; 0; 0; 0; 0; 0; x_sky(7)];
        case 'M2'
            xf_guess = [p.pc; p.r; 0; 0; 0; 0; x_sky(7)];
    end

    eta_guess = log(tf_guess);

    if strcmp(manifold,'M1x')
        param_guess = [eta_guess; 1; 1; 1; 1];   % eta + 4 multipliers
    else
        param_guess = [eta_guess; 1; 1; 1];      % eta + 3 multipliers
    end

    solinit = bvpinit(s_mesh, @(s) guess_fun(s,x_sky,xf_guess), param_guess);

    opts = bvpset('RelTol',1e-5,'AbsTol',1e-7,'NMax',5000);

    sol = bvp4c(@(s,y,param) phaseB_ode(s,y,param,p), ...
                @(ya,yb,param) phaseB_bc(ya,yb,param,p,x_sky,manifold), ...
                solinit,opts);

    s = linspace(0,1,400);
    y = deval(sol,s).';

    eta = sol.parameters(1);
    tf = exp(eta);
    t = s.'*tf;

    x = y(:,1:7);
    lam = y(:,8:14);

    u = zeros(length(t),2);
    H = zeros(length(t),1);

    for k = 1:length(t)
        [u(k,:),H(k)] = phaseB_control_H(x(k,:).',lam(k,:).',p);
    end

    Tphys = p.mh*p.g + u(:,1);

    integrand = zeros(length(t),1);
    for k = 1:length(t)
        uk = u(k,:).';
        integrand(k) = uk.'*p.R*uk;
    end

    J = trapz(t,integrand);

    fprintf('\nResults for %s:\n',manifold);
    fprintf('tf                    = %.6f s\n',tf);
    fprintf('J                     = %.6e\n',J);
    fprintf('H(tf)                 = %.6e\n',H(end));
    fprintf('Hamiltonian variation = %.6e\n',max(H)-min(H));
    fprintf('min physical T        = %.6f N\n',min(Tphys));
    fprintf('max physical T        = %.6f N\n',max(Tphys));
    fprintf('max abs tau           = %.6f N*m\n',max(abs(u(:,2))));
    fprintf('min altitude          = %.6f m\n',min(x(:,2)));

    if min(Tphys) < p.Tmin || max(Tphys) > p.Tmax
        fprintf('WARNING: thrust bounds violated.\n');
    else
        fprintf('Thrust bounds satisfied.\n');
    end

    if max(abs(u(:,2))) > p.taumax
        fprintf('WARNING: torque bounds violated.\n');
    else
        fprintf('Torque bounds satisfied.\n');
    end

    % Terminal manifold residual
    switch manifold
        case 'M1x'
            psi = [x(end,1)-p.land_px;
                   x(end,2);
                   x(end,4);
                   x(end,5)];
    
        case 'M2'
            psi = [(x(end,1)-p.pc)^2 + x(end,2)^2 - p.r^2;
                    x(end,4);
                    x(end,5)];
    end



    fprintf('terminal manifold residual = %.6e\n',norm(psi));

    res.manifold = manifold;
    res.t = t;
    res.x = x;
    res.lam = lam;
    res.u = u;
    res.H = H;
    res.Tphys = Tphys;
    res.J = J;
    res.tf = tf;

    plot_phaseB_result(res,p);
end

function yguess = guess_fun(s,x0,xf)

    xguess = (1-s)*x0 + s*xf;

    % Costate guess
    lamguess = [1; 1; 0.5; 1; 1; 0.5; 0];

    yguess = [xguess; lamguess];
end

function dyds = phaseB_ode(~,y,param,p)

    eta = param(1);
    tf = exp(eta);

    x = y(1:7);
    lam = y(8:14);

    [u,~] = phaseB_control_H(x,lam,p);

    xdot = p.A*x + p.B*u.';
    lamdot = -p.A.'*lam;

    dyds = tf*[xdot; lamdot];
end

function bc = phaseB_bc(ya,yb,param,p,x_sky,manifold)

    eta = param(1);
    tf = exp(eta);
    mu = param(2:end);

    xa = ya(1:7);
    xf = yb(1:7);
    lamf = yb(8:14);

    [~,Hf] = phaseB_control_H(xf,lamf,p);

    switch manifold
        case 'M1x'
            psi = [xf(1)-p.land_px; xf(2); xf(4); xf(5)];
        
            Dpsi = [1 0 0 0 0 0 0;
                    0 1 0 0 0 0 0;
                    0 0 0 1 0 0 0;
                    0 0 0 0 1 0 0];

        case 'M2'
            h1 = (xf(1)-p.pc)^2 + xf(2)^2 - p.r^2;
            psi = [h1; xf(4); xf(5)];

            Dpsi = [2*(xf(1)-p.pc), 2*xf(2), 0, 0, 0, 0, 0;
                    0, 0, 0, 1, 0, 0, 0;
                    0, 0, 0, 0, 1, 0, 0];
    end

    transversality = lamf - Dpsi.'*mu;

    bc = [
        xa - x_sky;
        psi;
        transversality;
        Hf
    ];
end

function [u,H] = phaseB_control_H(x,lam,p)

    ucol = -0.5*(p.R\(p.B.'*lam));

    % Unconstrained PMP control.
    % For physical validation, we check bounds afterward.
    u = ucol.';

    H = ucol.'*p.R*ucol + p.rho + lam.'*(p.A*x + p.B*ucol);
end

function plot_phaseB_result(res,p)

    figure('Name',['Phase B ' res.manifold ' States']);
    labels = {'p_x','p_z','v_x','v_z','theta','omega','m'};
    tiledlayout(4,2);

    for k = 1:7
        nexttile;
        plot(res.t,res.x(:,k),'LineWidth',1.5);
        grid on;
        xlabel('Time [s]');
        ylabel(labels{k});
        title([res.manifold ' state: ' labels{k}]);
    end

    figure('Name',['Phase B ' res.manifold ' Controls']);
    tiledlayout(3,1);

    nexttile;
    plot(res.t,res.Tphys,'LineWidth',1.5); hold on;
    yline(p.Tmin,'--');
    yline(p.Tmax,'--');
    grid on;
    ylabel('T [N]');
    title([res.manifold ' Physical Thrust']);

    nexttile;
    plot(res.t,res.u(:,2),'LineWidth',1.5); hold on;
    yline(p.taumin,'--');
    yline(p.taumax,'--');
    grid on;
    ylabel('\tau [N*m]');
    title([res.manifold ' Torque']);

    nexttile;
    plot(res.t,res.H,'LineWidth',1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('H');
    title([res.manifold ' Hamiltonian']);

    figure('Name',['Phase B ' res.manifold ' Trajectory with Thrust Vectors']);
    plot(res.x(:,1),res.x(:,2),'LineWidth',2); hold on; grid on; axis equal;
    xlabel('p_x [m]');
    ylabel('p_z [m]');
    title([res.manifold ' Landing Trajectory with Thrust Vectors']);
    
    if strcmp(res.manifold,'M2')
        ang = linspace(0,2*pi,200);
        xc = p.pc + p.r*cos(ang);
        zc = 0    + p.r*sin(ang);

        plot(xc,zc,'r--','LineWidth',2);

        % highlight landing platform
        patch(xc,zc,[1 0.85 0.85],'FaceAlpha',0.25,'EdgeColor','r','LineWidth',1.5);

        legend('Trajectory','Circular platform');
    end

    % Plot thrust direction arrows
    idx = round(linspace(1,length(res.t),12));

    for k = idx
        px = res.x(k,1);
        pz = res.x(k,2);
        th = res.x(k,5);
        Tmag = res.Tphys(k);

        % Thrust direction in inertial frame
        ux = sin(th);
        uz = cos(th);

        % Scale for plotting
        scale = 0.02*Tmag;

        quiver(px,pz,scale*ux,scale*uz,0,'LineWidth',1.2,'MaxHeadSize',1.5);
    end
end

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