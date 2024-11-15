clear; clc;

kp = 10;
kd = 2 * sqrt(kp) * 0.5;
alpha = 1;

f = @(x) [x(2); x(3); 0];
g = @(x) [0; 0; 1];
fz = @(z) 0;
gz = @(z) 1;
proj = @(x) x(1);

k = @(x, v) kp * (v - x(2)) - kd * x(3);

hd = @(z) z;

hd_filt = @(z, v_d) max(v_d, -alpha * hd(z));

%% Look at IC
x0 = [1, -1, -2];
vd = @(t) sign(sin(-2 * t)) - 0.5;
% vd = @(t) -1;

[t, x] = ode45(@(t, x) f(x) + g(x) * k(x,  hd_filt(proj(x), vd(t))), linspace(0, 100, 1000), x0);


figure(1)
clf
plot(t, x(:, 1))
yline(0)

figure(2)
clf
vdt = zeros(size(x, 1), 1);
filt_vdt = zeros(size(x, 1), 1);
for ii = 1:size(vdt, 1)
    vdt(ii) = vd(t(ii));
    filt_vdt(ii) = hd_filt(proj(x(ii, :)), vd(t(ii)));
end

subplot(2,1,1)
hold on
plot(t, vdt)
plot(t, filt_vdt)
plot(t, x(:, 2))
subplot(2,1,2)
plot(t, abs(filt_vdt - x(:, 2)))

%% Alright lets make some un-readable plots

kp = 10;
alpha = 1;

f = @(x) [x(2); 0];
g = @(x) [0; 1];
fz = @(z) 0;
gz = @(z) 1;
proj = @(x) x(1);

k = @(x, v) kp * (v - x(2));

hd = @(z) z;

hd_filt = @(z, v_d) max(v_d, -alpha * hd(z));

set(groot, 'DefaultAxesFontSize', 17); % Set default font size for axes labels and ticks 
set(groot, 'DefaultTextFontSize', 17); % Set default font size for text objects
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Set interpreter for axis tick labels 
set(groot, 'DefaultTextInterpreter', 'latex'); % Set interpreter for text objects (e.g., titles, labels) 
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultFigureRenderer', 'painters');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 15);
np = 10;
ndp = 11;
p_linspace = linspace(0, 4, np);
dp_linspace = linspace(-2, 2, ndp);

ps = zeros(np * ndp, 1);
dps = zeros(np * ndp, 1);
vds = zeros(np * ndp, 2);
vs = zeros(np * ndp, 2);
vd_const = -1;
fom_dyn = zeros(np * ndp, 2);
h_viol = zeros(np * ndp, 1);
v_err_signed = zeros(np * ndp, 1);

for i = 1:np
    for j = 1:ndp
        ps(np * (i-1) + j) = p_linspace(i);
        dps(np * (i-1) + j) = dp_linspace(j);
        x = [p_linspace(i); dp_linspace(j)];
        vds(np * (i-1) + j, :) = [vd_const; 0];
        vs(np * (i-1) + j, :) = [hd_filt(proj(x), vd_const); 0];
        fom_dyn(np * (i-1) + j, :) = f(x) + g(x) * k(x,  hd_filt(proj(x), vd_const));
    end
end

figure(1)
clf
hold on
quiver(ps, dps, vds(:, 1), vds(:, 2), 0.5)
quiver(ps, dps, vs(:, 1), vs(:, 2), 0.5)
quiver(ps, dps, fom_dyn(:, 1), fom_dyn(:, 2), 0.5)
xlabel('$p$')
ylabel('$\dot{p}$')
legend('$v_d$', '$v$', '$\nabla h f_{cl}(x)$')

% Define the grid for p and dp
p_linspace = linspace(0, 4, np * 10);
dp_linspace = linspace(-2, 2, ndp * 10);
[p_grid, dp_grid] = meshgrid(p_linspace, dp_linspace);

% Pre-allocate matrices for h_viol and v_err_signed
h_viol_mat = zeros(size(p_grid));
v_err_signed_mat = zeros(size(p_grid));

% Populate h_viol and v_err_signed on the grid
for i = 1:size(p_linspace, 2)
    for j = 1:size(dp_linspace, 2)
        x = [p_grid(j, i); dp_grid(j, i)];
        v = hd_filt(proj(x), vd_const);
        fom_dyn = f(x) + g(x) * k(x, v);
        
        % Calculate the values at each grid point
        h_viol_mat(j, i) = -min(0, fom_dyn(1) + alpha * p_grid(j, i));
        v_err_signed_mat(j, i) = v - fom_dyn(1);
    end
end

% Plot using surf
figure(2)
clf
subplot(1,2,1)
surf(p_grid, dp_grid, h_viol_mat, FaceColor='interp')
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("CBF Violation")
colorbar
view(0, 90);

subplot(1,2,2)
surf(p_grid, dp_grid, abs(v_err_signed_mat), FaceColor='interp', EdgeColor='none')
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("Signed Velocity Tracking Error")
colorbar
view(0, 90);

figure(3)
clf
mu = 0.8;
subplot(1,2,1)
contour(p_grid, dp_grid, p_grid - mu * v_err_signed_mat.^2, [0 0])
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("Candidate BarF")
% colorbar
subplot(1,2,2)
surf(p_grid, dp_grid, p_grid - mu * v_err_signed_mat.^2, FaceColor='interp', EdgeColor='none')
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("Candidate BarF")
colorbar

figure(4)
clf
mu = 0.8;
subplot(1,2,1)
contour(p_grid, dp_grid, p_grid - mu * h_viol_mat, [0 0])
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("Candidate BarF")
% colorbar
subplot(1,2,2)
surf(p_grid, dp_grid, p_grid - mu * h_viol_mat, FaceColor='interp', EdgeColor='none')
xlabel('$p$')
ylabel('$\dot{p}$')
zlabel("Candidate BarF")
colorbar
