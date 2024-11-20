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
kp = 5;
v_max = 2;
function v = v_desired(z, v_max, kp)
v = -kp * z;
if norm(v) > v_max
    v = v / norm(v) * v_max;
end
end
v_des = @(z) v_desired(z, v_max, kp);

kd = 2;
k_track = @(x, v) -kd * [x(3) - v(1); x(4) - v(2)];

fx = @(x, u) [x(3); x(4); u(1); u(2)];
fz = @(z) zeros(2, 1);
gz = @(z) eye(2);

c1 = [-0.5; -2];
c2 = [-2; -3];
r = 0.5;
rho = 2;
% h = @(z) -1 / rho * log(exp(-rho * (norm(z - c1) - r)) + exp(-rho * (norm(z - c2) - r)));
% Jh = @(z) -1 / rho / (exp(-rho * (norm(z - c1) - r)) + exp(-rho * (norm(z - c2) - r))) * ...
%     (-rho * (z - c1)' / norm(z - c1) * exp(-rho * (norm(z - c1) - r)) -rho * (z - c2)' / norm(z - c2) * exp(-rho * (norm(z - c2) - r)));
h = @(z) min(norm(z - c1) - r, norm(z - c2) - r);

function jh_ = Jh_func(z, c1, c2, r)
if norm(z - c1) - r < norm(z - c2) - r
    jh_ = (z - c1)' / norm(z - c1);
else
    jh_ = (z - c2)' / norm(z - c2);
end
end
Jh = @(z) Jh_func(z, c1, c2, r);

z1 = linspace(-4, 0);
z2 = linspace(-5, 0);

[Z1, Z2] = meshgrid(z1, z2);
hZ = zeros(size(Z1));
for i = 1:numel(Z1)
    hZ(i) = h([Z1(i); Z2(i)]);
end
% figure(3)
% clf
% hold on
% surf(Z1, Z2, hZ, 'FaceColor','interp', 'EdgeColor','none')
% % contour(Z1, Z2, hZ,[0 0])
% colorbar
% theta = linspace(0, 2 * pi);
% plot3(c1(1) + r * sin(theta), c1(2) + r * cos(theta), zeros(size(theta)), 'r')
% plot3(c2(1) + r * sin(theta), c2(2) + r * cos(theta), zeros(size(theta)), 'r')
% axis square
% axis equal
% xlabel('X')
% ylabel('Y')
% zlabel('h')
% view([0 90])

% figure(4)
% clf
% hold on
% plot(c1(1) + r * sin(theta), c1(2) + r * cos(theta), 'r')
% plot(c2(1) + r * sin(theta), c2(2) + r * cos(theta), 'r')
% x1s = [0, 0, 0, -1, -1, -1, -2, -2, -2, -3, -3, -3];
% y1s = [-0.5, -2.5, -4, -1, -2.5, -4, -1, -2.2, -3.8, -2, -3, -4];
% Jh1s = [];
% Jh2s = [];
% for ii = 1:numel(x1s)
%     Jh1 = Jh([x1s(ii); y1s(ii)]);
%     Jh1s = [Jh1s, Jh1(1)];
%     Jh2s = [Jh2s, Jh1(2)];
% end
% quiver(x1s, y1s, Jh1s, Jh2s)
% axis square
% axis equal

Pi = @(x) [x(1); x(2)];
JPi = @(x) [1 0 0 0; 0 1 0 0];

alpha = 5;
epsilon = 0;
T = 4;
dt = 0.01;
decimation = 5;
iters = 15;
tol = 1e-5;
K_delta = 1;
verbose = 1;

smooth_psf = smoother_predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, JPi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);
psf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);
sf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, false);

% x0 = [-3, -5, 0, 0];
x0 = [-2.1, -4, 0, 2];

%% Default Experiment
tic;
[t2, x2, h2, vd2, v2, d2] = smooth_psf.apply_cl(x0, 10, 0.01);
[t, x, h, vd, v, d] = psf.apply_cl(x0, 10, 0.01);
[t1, x1, h1, vd1, v1, d1] = sf.apply_cl(x0, 10, 0.01);

fprintf("Runtime: %0.2f\n", toc)

%% Plot to visualize
figure(1)
clf
subplot(2,2,1)

hold on
plot(t, h1)
plot(t, h)
% if any(h < 0)
%     plot(t(h < 0), h(h < 0), 'ro')
% end

% if any(h1 < 0)
%     plot(t1(h1 < 0), h1(h1 < 0), 'go')
% end
plot(t, h2)
xlabel('Time (s)')
ylabel('$h(\Pi(x))$')
legend( 'Nominal SF', 'Predictive SF','Barrier PSF', AutoUpdate=false)
yline(0, 'k')
% sgtitle(sprintf('T: %d, Iters: %d, alpha: %d, epsilon: %d, Kp: %0.2f, Kd: %0.2f     Max Violation: %0.2e', T, iters, alpha, epsilon, kp, kd, -min(0, min(h))))
sgtitle(sprintf('Maximum Violations: Nominal: %0.2e PSF %0.2e BPSF: %0.2e', -min(0, min(h1)), -min(0, min(h)), -min(0, min(h2))))
subplot(2,2,2)
hold on
plot(x1(:, 1), x1(:, 2))
plot(x(:, 1), x(:, 2))
plot(x2(:, 1), x2(:, 2))
legend('Nominal SF', 'Predictive SF', 'Barrier PSF', AutoUpdate=false)
plot(x0(1), x0(2), 'rx')
plot(0, 0, 'go')
axis equal
theta = linspace(0, 2 * pi);
plot(c1(1) + r * sin(theta), c1(2) + r * cos(theta), 'r')
plot(c2(1) + r * sin(theta), c2(2) + r * cos(theta), 'r')

xlabel('X (m)')
ylabel('Y (m)')

% subplot(2,2,3)
% hold on
% plot(t, vd(:, 1), LineWidth=4)
% plot(t, v(:, 1))
% plot(t, x(:, 3))
% legend('$v_x$ desired', '$v_x$ filtered', '$v_x$ actual')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% 
% subplot(2,2,4)
% hold on
% plot(t, vd(:, 2))
% plot(t, v(:, 2))
% plot(t, x(:, 4))
% legend('$v_z$ desired', '$v_z$ filtered', '$v_z$ actual')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')

subplot(2,2,3)
hold on
plot(t, vd(:, 1), LineWidth=4)
plot(t, v(:, 1))
plot(t, x(:, 3))
legend('$v_x$ desired', '$v_x$ filtered', '$v_x$ actual')
xlabel('Time (s)')
ylabel('Velocity (m/s)')

subplot(2,2,4)
hold on
plot(t, d1)
plot(t, d)
plot(t, d2)
legend('Nominal SF', 'Predictive SF', 'Barrier PSF', AutoUpdate=false)
xlabel('Time (s)')
ylabel('Delta')