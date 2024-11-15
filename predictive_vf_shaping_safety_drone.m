clear; clc;

% Define tracking controller
kp = 10;
% kd = 2;
kd = 2 * sqrt(10);
kv = 2;
kp_rom = 1;


% Define planner/tracker dynamics
function [f, tau] = vel_controller(x, v, m, I, g, kv, kp, kd)
a_des = kv * (v - x(4:5));
f_des = m * (a_des + [0; g]);
f = f_des' * [sin(x(3)); cos(x(3))];
theta_des = atan2(f_des(1), f_des(2));
tau = kp * (theta_des - x(3)) - kd * x(6);
end
function dotx = fx(x, v, k, m, I, g)
[f, tau] = k(x, v);
dotx = [x(4); x(5); x(6); f / m * sin(x(3)); -g + f / m * cos(x(3)); tau / I];
end

m = 1;
I = 1;
g = 9.81;
k = @(x_, v_) vel_controller(x_, v_, m, I, g, kv, kp, kd);
%% Test vel controller
x0 = [0; 0; 0; 0; 0; 0];
vd = [-1; -2];

[t, x] = ode45(@(t, x) fx(x, vd, k, m, I, g), [0, 10], x0);

figure(1)
clf
subplot(2,1,1)
plot(t, x(:, 1:3))
legend('x', 'y', 'theta')
xlabel('Time')
ylabel('State')
subplot(2,1,2)
plot(t, x(:, 4:6))
hold on
plot(t, ones(size(t)) * vd')
legend('dx', 'dy', 'dtheta', 'vdx', 'vdy')
xlabel('Time')
ylabel('State')


%% Safety stuff
global delta

% Predictive safety filter
function v = predictive_safety_filter(t, vd, x0, k, alpha, epsilon, T, iters, m, I, g)
    global delta
    fprintf("t: %0.4f\n", t)
    fprintf("\tI: %d, delta: %0.2f\n", 0, delta)
    for i = 1:iters
        % Roll out planner tracker with buffered safety function
        h_d = @(x) x(2) - delta;
        rom_filt = @(vd, x) max(vd(x(1:2)), -alpha * h_d(x) + epsilon);
        [t, x] = ode45(@(t, x) fx(x, rom_filt(vd, x), k, m, I, g), [0, T], x0);

        % Evaluate barrier function on full order model
        h_x = x(:, 2);
        h_bar = min(h_x);

        % Buffer barrier function by violation of FoM
        delta = delta - h_bar;
        fprintf("\tI: %d, delta: %0.2f\n", i, delta)
        % set(gca,'ColorOrderIndex',1)
        % plot(t, y(:, 3))
        % plot(t, y(:, 1))
        % yline(0, 'k')
        % drawnow
    end  % Hope for convergence
    h_d = @(x) x(2) - max(0, delta);
    rom_filt = @(vd, x) max(vd(x(1:2)), -alpha * h_d(x) + epsilon);
    v = rom_filt(vd, x0);
    % pause
end

%% Alright lets try it
T = 6;          % Horizon to integrate forwards
alpha = 1;     % class K function for h_dot
delta = 0;     % Initial guess for backoff term
x0 = [0; 2; pi/4; 0; 0; 0];
iters = 20;      % Iterations of predictive filter to run
vd = @(z) [-kp_rom * z(1); -5];        % Desired velocity for RoM
epsilon = 0;  % robustness term

tic;
[t, x] = ode45(@(t, x) fx(x, predictive_safety_filter(t, vd, x, k, alpha, epsilon, T, iters, m, I, g), k, m, I, g), [0, 10], x0);
fprintf("Runtime: %0.2f", toc)

% Plot to visualize
figure(1)
clf
plot(t, x(:, 2))
legend('y', AutoUpdate=false)
yline(0, 'k')
if any(x(:, 2) < 0)
    xline(t(x(:, 2) < 0), 'r', LineWidth=2, Alpha=0.5)
end
xlabel('Time (s)')
ylabel('Vertical Position (m)')
title(sprintf('T: %d, Iters: %d, alpha: %d, Kp: %0.2f, Kd: %0.2f', T, iters, alpha, kp, kd))

figure(2)
clf
hold on
plot(x(:, 1), x(:, 2))
yline(0, 'k')
if any(x(:, 2) < 0)
    scatter(x(x(:, 2) < 0, 1), x(x(:, 2) < 0, 2), 'ro')
end
xlabel('X (m)')
ylabel('Z (m)')
title(sprintf('T: %d, Iters: %d, alpha: %d, Kp: %0.2f, Kd: %0.2f', T, iters, alpha, kp, kd))