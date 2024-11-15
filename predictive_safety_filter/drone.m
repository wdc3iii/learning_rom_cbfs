clear; clc;
set(groot, 'DefaultAxesFontSize', 17); % Set default font size for axes labels and ticks
set(groot, 'DefaultTextFontSize', 17); % Set default font size for text objects
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Set interpreter for axis tick labels
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultFigureRenderer', 'painters');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 15);

% Define all the stuff
kp_rom = 1;
kp = 20;
kd = 2 * sqrt(10);
kv = 2;

% Define planner/tracker dynamics
function u = vel_controller(x, v, m, I, g, kv, kp, kd)
a_des = kv * (v - [x(4); x(5)]);
f_des = m * (a_des + [0; g]);
f = f_des' * [sin(x(3)); cos(x(3))];
theta_des = atan2(f_des(1), f_des(2));
tau = kp * (theta_des - x(3)) - kd * x(6);
u = [f tau];
end

function dotx = f_x(x, u, m, I, g)
dotx = [x(4); x(5); x(6); u(1) / m * sin(x(3)); -g + u(1) / m * cos(x(3)); u(2) / I];
end

m = 1;
I = 1;
g = 9.81;
v_des = @(z) [-kp_rom * z(1); -5];
k_track = @(x_, v_) vel_controller(x_, v_, m, I, g, kv, kp, kd);

fx = @(x, u) f_x(x, u, m, I, g);
fz = @(z) zeros(2, 1);
gz = @(z) eye(2);
h = @(z) z(2);
Jh = @(z) [0 1];
Pi = @(x) [x(1); x(2)];

alpha = 1;
epsilon = 0;
T = 6;
dt = 0.01;
decimation = 5;
iters = 15;
tol = 1e-5;
K_delta = 0.5;
verbose = 1;

psf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);

x0 = [0 1 pi/4 0 0 0];

%% Default Experiment
tic;
[t, x, h, vd, v, d] = psf.apply_cl(x0, 20, 0.01);
fprintf("Runtime: %0.2f\n", toc)

%% Plot to visualize
figure(1)
clf
subplot(2,2,1)
hold on
plot(t, h)
yline(0, 'k')
if any(h < 0)
    plot(t(h < 0), h(h < 0), 'ro')
end
xlabel('Time (s)')
ylabel('$h(\Pi(x))$')

subplot(2,2,2)
hold on
plot(x(:, 1), x(:, 2))
yline(0, 'k')
if any(x(:, 2) < 0)
    plot(x(x(:, 2) < 0, 1), x(x(:, 2) < 0, 2), 'ro')
end
xlabel('X (m)')
ylabel('Z (m)')
sgtitle(sprintf('T: %d, Iters: %d, alpha: %d, epsilon: %d, Kp: %0.2f, Kd: %0.2f     Max Violation: %0.2e', T, iters, alpha, epsilon, kp, kd, -min(0, min(x(:, 2)))))

subplot(2,2,3)
hold on
plot(t, vd(:, 1), LineWidth=4)
plot(t, v(:, 1))
plot(t, x(:, 4))
legend('$v_x$ desired', '$v_x$ filtered', '$v_x$ actual')

subplot(2,2,4)
hold on
plot(t, vd(:, 2))
plot(t, v(:, 2))
plot(t, x(:, 5))
legend('$v_z$ desired', '$v_z$ filtered', '$v_z$ actual')

% %% Max Violation as a function of T, iters
% 
% Ts = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5];
% iters = [1 2 3 4 5];
% 
% [Ts, Is] = meshgrid(Ts, iters);
% 
% maxViols = zeros(size(Ts));
% 
% parfor i = 1:numel(Ts)
%     psf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, Ts(i), dt, Is(i), tol, K_delta, verbose, 2);
%     [t, x] = psf.apply_cl(x0, 10, 0.01);
% 
%     % Plot to visualize
%     figure(1)
%     clf
%     hold on
%     plot(t, x(:, 1))
%     yline(0, 'k')
%     if any(x(:, 1) < 0)
%         plot(t(x(:, 1) < 0), x(x(:, 1) < 0, 1), 'ro')
%     end
%     xlabel('Time (s)')
%     ylabel('$h(\Pi(x))$')
%     max_viol = -min(0, min(x(:, 1)));
%     maxViols(i) = max_viol;
%     title(sprintf('T: %0.2f, Iters: %d, alpha: %d, epsilon: %d, Kp: %0.2f, Kd: %0.2f     Max Violation: %0.2e', Ts(i), Is(i), alpha, epsilon, kp, kd, max_viol))
%     drawnow
% end
% 
% figure(2)
% clf
% surf(Ts, Is, log(maxViols), 'FaceColor','interp')
% colorbar
% xlabel('Horizon (s)')
% ylabel('Iterations')
% view([0 90])
% title('Log Max Violation')