clear; clc;
set(groot, 'DefaultAxesFontSize', 17); % Set default font size for axes labels and ticks
set(groot, 'DefaultTextFontSize', 17); % Set default font size for text objects
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Set interpreter for axis tick labels
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultFigureRenderer', 'painters');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 15);

% Define nominal controller
kp = 1;
v_max = 2;
function v = v_desired(z, v_max, kp)
v = -kp * z;
if norm(v) > v_max
    v = v / norm(v) * v_max;
end
end
v_des = @(z) v_desired(z, v_max, kp);

% Define tracking controller
kd = 1;
k_track = @(x, v) -kd * [x(3) - v(1); x(4) - v(2)];

% Define dynamics
fx = @(x, u) [x(3); x(4); u(1); u(2)];
fz = @(z) zeros(2, 1);
gz = @(z) eye(2);

% Define safe set
c1 = [-0.5; -2];
c2 = [-2; -3];
r = 0.5;
h = @(z) min(norm(z - c1) - r, norm(z - c2) - r);
function jh_ = Jh_func(z, c1, c2, r)
if norm(z - c1) - r < norm(z - c2) - r
    jh_ = (z - c1)' / norm(z - c1);
else
    jh_ = (z - c2)' / norm(z - c2);
end
end
Jh = @(z) Jh_func(z, c1, c2, r);

% Define projection
Pi = @(x) [x(1); x(2)];
JPi = @(x) [1 0 0 0; 0 1 0 0];

% Define algorithm parameters
alpha = 5;
epsilon = 0;
T = 5;
dt = 0.01;
iters = 20;
tol = 1e-5;
K_delta = 1;
verbose = 1;

sf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, false);
bpsf = smoother_predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, JPi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);

x0 = [-3, -5, 0, 0];

%% Default Experiment
[tsf, xsf, hsf, vdsf, vsf, dsf] = sf.apply_cl(x0, 10, 0.01);
[tb, xb, hb, vdb, vb, db] = bpsf.apply_cl(x0, 10, 0.01);

%% Plot to visualize

load("/home/wcompton/repos/legged_gym_dev/predictive_cbfs/predictive_cbfs/double_single_int_j4fo3cmm/eval_data.mat")
t = [0:size(delta, 2) - 1] *  0.01;
figure(1)
clf
subplot(1,2,1)

hold on
plot(tsf, dsf)
plot(tb, db)
plot(t, delta)
xlabel('Time (s)')
ylabel('$\delta$')
legend('Nominal SF', 'Predictive SF', 'Learned PSF', AutoUpdate=false)
fprintf('Maximum Violations: Nominal: %0.2e PSF: %0.2e, Learned PSF: %0.2e', -min(0, min(hsf)), -min(0, min(hb)), -min(0, min(h)))

subplot(1,2,2)
hold on
plot(x0(1), x0(2), 'rx')
plot(0, 0, 'go')
axis equal
theta = linspace(0, 2 * pi);
rectangle('Position', [c1(1) - r, c1(2) - r, 2*r, 2*r], 'Curvature', [1, 1], 'FaceColor', 'k', 'EdgeColor', 'None');
rectangle('Position', [c2(1) - r, c2(2) - r, 2*r, 2*r], 'Curvature', [1, 1], 'FaceColor', 'k', 'EdgeColor', 'None');
set(gca,'ColorOrderIndex',1)
plot(xsf(:, 1), xsf(:, 2), DisplayName='Nominal SF')
plot(xb(:, 1), xb(:, 2), DisplayName='Predictive SF')
plot(x(:, 1), x(:, 2), DisplayName='Learned PSF')
plts = findobj(gca, 'Type', 'line');
legend(plts([3 2 1]));
xlabel('X (m)')
ylabel('Y (m)')