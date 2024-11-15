clear; clc;

% Define tracking controller
kp = 10;
kd = 2;
% kd = 2 * sqrt(10);
k = @(x, v) -kp * (x(2) - v) - kd * x(3);

% Define planner/tracker dynamics
function dotx = fx(x, v, k)
dotx = [x(2); x(3); k(x, v)];
end

global delta

% Predictive safety filter
function v = predictive_safety_filter(t, vd, x0, k, alpha, epsilon, T, iters)
    global delta
    figure(1)
    clf
    hold on
    fprintf("t: %0.4f\n", t)
    fprintf("\tI: %d, delta: %0.2f\n", 0, delta)
    for i = 1:iters
        % Roll out planner tracker with buffered safety function
        h_d = @(x) x(1) - delta;
        rom_filt = @(vd, x) max(vd, -alpha * h_d(x) + epsilon);
        [t, x] = ode45(@(t, x) fx(x, rom_filt(vd, x), k), [0, T], x0);

        % Evaluate barrier function on full order model
        h_x = x(:, 1);
        h_bar = min(h_x);

        % Buffer barrier function by violation of FoM
        delta = delta - h_bar;
        % v = rom_filt(vd, x0);
        fprintf("\tI: %d, delta: %0.2f\n", i, delta)
        % set(gca,'ColorOrderIndex',1)
        % plot(t, y(:, 3))
        % plot(t, y(:, 1))
        % yline(0, 'k')
        % drawnow
    end  % Hope for convergence
    h_d = @(x) x(1) - max(0, delta);
    rom_filt = @(vd, x) max(vd, -alpha * h_d(x) + epsilon);
    v = rom_filt(vd, x0);
    % pause
end

%% Alright lets try it
T = 6;          % Horizon to integrate forwards
alpha  = 1;     % class K function for h_dot
delta = 0;     % Initial guess for backoff term
x0 = [1; 0; 0];    % FoM IC
z0 = 1;         % RoM IC
iters = 20;      % Iterations of predictive filter to run
vd = -5;        % Desired velocity for RoM
epsilon = 0.1;  % robustness term

tic;
[t, x] = ode45(@(t, x) fx(x, predictive_safety_filter(t, vd, x, k, alpha, epsilon, T, iters), k), [0, 20], x0);
fprintf("Runtime: %0.2f", toc)

% Plot to visualize
figure(1)
clf
plot(t, x(:, 1))
legend('x1', AutoUpdate=false)
yline(0, 'k')
if any(x(:, 1) < 0)
    xline(t(x(:, 1) < 0), 'r', LineWidth=2, Alpha=0.5)
end
xlabel('Time (s)')
ylabel('Position (m)')
title(sprintf('T: %d, Iters: %d, alpha: %d, Kp: %0.2f, Kd: %0.2f', T, iters, alpha, kp, kd))