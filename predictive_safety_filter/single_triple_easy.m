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
v_des = @(z) -5;
kp = 10;
kd = 2;
k_track = @(x, v) -kp * (x(2) - v) - kd * x(3);

fx = @(x, u) [x(2); x(3); u];
fz = @(z) 0;
gz = @(z) 1;
h = @(z) z;
Jh = @(z) 1;
Pi = @(x) x(1);
JPi = @(x) [1 0 0];

alpha = 1;
epsilon = 0;
T = 4;
dt = 0.01;
decimation = 5;
iters = 15;
tol = 1e-5;
K_delta = 0.5;
verbose = 1;

sf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, false);
psf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);
bpsf = smoother_predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, JPi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, 2, true);

x0 = [1 0 0];

%% Default Experiment
tic;
[tsf, xsf, hsf, vdsf, vsf, dsf] = sf.apply_cl(x0, 10, 0.01);
[t, x, h, vd, v, d] = psf.apply_cl(x0, 10, 0.01);
[tb, xb, hb, vdb, vb, db] = bpsf.apply_cl(x0, 10, 0.01);

fprintf("Runtime: %0.2f\n", toc)

%% Plot to visualize
figure(1)
clf
subplot(1,2,1)
hold on
plot(tsf, hsf)
plot(t, h)
plot(tb, hb)
legend( 'Nominal SF', 'Predictive SF','Barrier PSF', AutoUpdate=false)
yline(0, 'k')
% sgtitle(sprintf('T: %d, Iters: %d, alpha: %d, epsilon: %d, Kp: %0.2f, Kd: %0.2f     Max Violation: %0.2e', T, iters, alpha, epsilon, kp, kd, -min(0, min(h))))
sgtitle(sprintf('Maximum Violations: Nominal: %0.2e PSF %0.2e BPSF: %0.2e', -min(0, min(hsf)), -min(0, min(h)), -min(0, min(hb))))
% if any(h < 0)
%     plot(t(h < 0), h(h < 0), 'ro')
% end
xlabel('Time (s)')
ylabel('$h(\Pi(x))$')

subplot(1,2,2)
hold on
plot(tsf, dsf)
plot(t, d)
plot(tb, db)
legend( 'Nominal SF', 'Predictive SF','Barrier PSF', AutoUpdate=false)
xlabel('Time (s)')
ylabel('$\delta$')

%% Max Violation as a function of T, iters

% Ts = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5];
% iters = [1 2 3 4 5];
% 
% [Ts, Is] = meshgrid(Ts, iters);
% 
% maxViols = zeros(size(Ts));
% 
% parfor i = 1:numel(Ts)
%     psf = predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, h, Jh, alpha, epsilon, Ts(i), dt, Is(i), tol, K_delta, verbose, 2);
%     [t, x, h, vd, v] = psf.apply_cl(x0, 10, 0.01);
% 
%     % Plot to visualize
%     figure(1)
%     clf
%     hold on
%     plot(t, h)
%     yline(0, 'k')
%     if any(h < 0)
%         plot(t(h < 0), h(h < 0), 'ro')
%     end
%     xlabel('Time (s)')
%     ylabel('$h(\Pi(x))$')
%     max_viol = -min(0, min(h));
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