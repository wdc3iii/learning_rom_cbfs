clear; clc;

% Define tracking controller
% kp = 1;
% kd = .5;
kp = 1;
kd = 2 * sqrt(kp);
k = @(x, z, v) -kp * (x(1) - z) - kd * (x(2) - v);

% Define planner/tracker dynamics
function dotxz = fxz(x, z, v, k)
dotz = v;
dotx = [x(2); k(x, z, v)];
dotxz = [dotx; dotz];
end

global delta

% Predictive safety filter
function v = predictive_safety_filter(t, vd, x0, z0, k, alpha, T, iters)
    global delta
    figure(1)
    clf
    hold on
    fprintf("t: %0.4f\n", t)
    fprintf("\tI: %d, delta: %0.2f\n", 0, delta)
    for i = 1:iters
        % Roll out planner tracker with buffered safety function
        h_d = @(z) z - delta;
        rom_filt = @(vd, z) max(vd, -alpha * h_d(z));
        [t, y] = ode45(@(t, y) fxz(y(1:2), y(3), rom_filt(vd, y(3)), k), [0, T], [x0; z0]);

        % Evaluate barrier function on full order model
        h_x = y(:, 1);
        h_bar = min(h_x);

        % Buffer barrier function by violation of FoM
        delta = delta - h_bar;
        v = rom_filt(vd, z0);
        fprintf("\tI: %d, delta: %0.2f\n", i, delta)
        % set(gca,'ColorOrderIndex',1)
        % plot(t, y(:, 3))
        % plot(t, y(:, 1))
        % yline(0, 'k')
        % drawnow
    end  % Hope for convergence
    h_d = @(z) z - max(0, delta);
    rom_filt = @(vd, z) max(vd, -alpha * h_d(z));
    v = rom_filt(vd, z0);
    % pause
end

%% Alright lets try it
T = 1;          % Horizon to integrate forwards
alpha  = 1;     % class K function for h_dot
delta = 0;     % Initial guess for backoff term
x0 = [1; 0];    % FoM IC
z0 = 1;         % RoM IC
iters = 5;      % Iterations of predictive filter to run
vd = -5;        % Desired velocity for RoM

tic;
[t, y] = ode45(@(t, y) fxz(y(1:2), y(3), predictive_safety_filter(t, vd, y(1:2), y(3), k, alpha, T, iters), k), linspace(0, 10, 1000), [x0; z0]);
fprintf("Runtime: %0.2f", toc)
xt = y(:, 1:2);
zt = y(:, 3);

% Plot to visualize
figure(1)
clf
hold on
plot(t, zt)
plot(t, xt(:, 1))
legend('z', 'x1', AutoUpdate=false)
yline(0, 'k')
if any(xt(:, 1) < 0)
    xline(t(xt(:, 1) < 0), 'r', LineWidth=2, Alpha=0.5)
end
xlabel('Time (s)')
ylabel('Position (m)')
title(sprintf('T: %d, Iters: %d, alpha: %d, Kp: %0.2f, Kd: %0.2f', T, iters, alpha, kp, kd))

figure(2)
clf
hold on
plot(t, -alpha * xt(:, 1))
plot(t, xt(:, 2))
set(groot, 'DefaultLegendInterpreter', 'latex');
legend('$-\alpha h(\Pi(x))$', '$\dot{h}(x)$')
