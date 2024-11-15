clear;clc;
set(groot, 'DefaultAxesFontSize', 17); % Set default font size for axes labels and ticks
set(groot, 'DefaultTextFontSize', 17); % Set default font size for text objects
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex'); % Set interpreter for axis tick labels
set(groot, 'DefaultTextInterpreter', 'latex'); % Set interpreter for text objects (e.g., titles, labels)
set(groot, 'DefaultLegendInterpreter', 'latex');
set(groot, 'DefaultFigureRenderer', 'painters');
set(groot, 'DefaultLineLineWidth', 2);
set(groot, 'DefaultLineMarkerSize', 15);

vd = -5;
x0 = 1;
v0 = 0;
alpha = 1;
epsilon = 5;
% kp = 10;
% kd = 2 * sqrt(kp);
kp = 1;
kd = .5;
k = @(x, z, v) -kp * (x(1) - z) - kd * (x(2) - v);

fh = 1;

function dotxz = fxz(x, z, v, k)
dotz = v;
dotx = [x(2); k(x, z, v)];
dotxz = [dotx; dotz];
end


figure(fh)
clf
figure(fh + 1)
clf
sigmas = [0 1 2 3];
for ii = 1:size(sigmas, 2)

    sigma = sigmas(ii);
    H = @(z,x) z - sigma / 2 * (z - x(1)).^2;
    dotH = @(z, x, v) v - sigma * (z - x(1)) * (v - x(2));
    safety_filt = @(v_d, x, z) max(v_d, (-alpha * (z - sigma / 2 * (z - x(1)).^2) - sigma * (z - x(1)) * x(2) + 1/epsilon) / (1 - sigma * (z - x(1))));

    run_ic = @(x0, z0) ode45(@(t, y) fxz(y(1:2), y(3), safety_filt(vd, y(1:2), y(3)), k), [0, 10], [x0; z0]);

    [t, y] = run_ic([x0; v0], x0);
    xt = y(:, 1:2);
    zt = y(:, 3);

    figure(fh)
    subplot(size(sigmas, 2), 2, 2 * (ii - 1) + 1)
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
    title(sprintf('alpha: %d, sigma: %d, Kp: %d, Kd: %d', alpha, sigma, kp, kd))

    subplot(size(sigmas, 2), 2, 2 * ii)
    Ht = zeros(size(t));
    for tt = 1:size(t, 1)
        Ht(tt) = H(zt(tt), xt(tt, :));
    end
    plot(t, Ht)
    xlabel('Time (s)')
    ylabel('H(x, z)')

    figure(fh + 1)
    hold on
    vt = zeros(size(t));
    for tt = 1:size(t, 1)
        vt(tt) = safety_filt(vd, xt(tt, :), zt(tt));
    end
    plot(t, vt)
    ylim([-abs(vd), abs(vd)])
    drawnow
end


%% More noise
epsilon = 0;

figure(fh + 2)
clf
figure(fh + 3)
clf
sigmas = [0.1 1 2 3];
function v = safety_filt_w_dote(v_d, x, z, s, a)
    if -s * (z-x(1)) * v_d + (1 + s * (z - x(1))) * x(2) >= - a * (z - s * (z - x(1)).^2)
        v = v_d;
    else
        v = (-a * (z - s / 2 * (z - x(1)).^2) - (1 + s * (z - x(1)))) / (- s * (z - x(1)));
    end
end
for ii = 1:size(sigmas, 2)

    sigma = sigmas(ii);
    H = @(z,x) z - sigma / 2 * (z - x(1)).^2;
    dotH = @(z, x, v) v - sigma * (z - x(1)) * (v - x(2));
    
    run_ic = @(x0, z0) ode45(@(t, y) fxz(y(1:2), y(3), safety_filt_w_dote(vd, y(1:2), y(3), sigma, alpha), k), [0, 10], [x0; z0 - 0.5]);

    [t, y] = run_ic([x0; v0], x0);
    xt = y(:, 1:2);
    zt = y(:, 3);

    figure(fh + 2)
    subplot(size(sigmas, 2), 2, 2 * (ii - 1) + 1)
    hold on
    plot(t, zt)
    plot(t, xt(:, 1))
    legend('z', 'x1', AutoUpdate=false)
    yline(0, 'k')
    % if any(xt(:, 1) < 0)
    %     xline(t(xt(:, 1) < 0), 'r', LineWidth=2, Alpha=0.5)
    % end
    xlabel('Time (s)')
    ylabel('Position (m)')
    title(sprintf('alpha: %d, sigma: %d, Kp: %d, Kd: %d', alpha, sigma, kp, kd))

    subplot(size(sigmas, 2), 2, 2 * ii)
    Ht = zeros(size(t));
    for tt = 1:size(t, 1)
        Ht(tt) = H(zt(tt), xt(tt, :));
    end
    plot(t, Ht)
    xlabel('Time (s)')
    ylabel('H(x, z)')

    figure(fh + 3)
    hold on
    vt = zeros(size(t));
    for tt = 1:size(t, 1)
        vt(tt) = safety_filt_w_dote(vd, xt(tt, :), zt(tt), sigma, alpha);
    end
    plot(t, vt)
    ylim([-abs(vd), abs(vd)])
    drawnow
end

%% More noise
epsilon = 0;
delt = 0.9;

figure(fh + 4)
clf
figure(fh + 5)
clf
sigmas = [0.1 1 2 3];
function v = safety_filt_w_norm_dote(v_d, x, z, s, a, delt)
    term = -a * (z - s * (z - x(1)).^2) - s * (z - x(1)) * x(2);
    v1 = (term - delt * x(2)) / (1 - delt -s * (z - x(1)));
    v2 = (term + delt * x(2)) / (1 + delt -s * (z - x(1)));
    if v1 - x(2) >= 0
        v = max(v_d, v1);
    elseif v2 - x(2) <= 0
        v = max(v_d,  v2);
    else
        v = nan;
    end
end
for ii = 1:size(sigmas, 2)

    sigma = sigmas(ii);
    H = @(z,x) z - sigma / 2 * (z - x(1)).^2;
    dotH = @(z, x, v) v - sigma * (z - x(1)) * (v - x(2));
    
    run_ic = @(x0, z0) ode45(@(t, y) fxz(y(1:2), y(3), safety_filt_w_norm_dote(vd, y(1:2), y(3), sigma, alpha, delt), k), [0, 10], [x0; z0 - 0.1]);

    [t, y] = run_ic([x0; v0], x0);
    xt = y(:, 1:2);
    zt = y(:, 3);

    figure(fh + 4)
    subplot(size(sigmas, 2), 2, 2 * (ii - 1) + 1)
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
    title(sprintf('alpha: %d, sigma: %d, Kp: %d, Kd: %d', alpha, sigma, kp, kd))

    subplot(size(sigmas, 2), 2, 2 * ii)
    Ht = zeros(size(t));
    for tt = 1:size(t, 1)
        Ht(tt) = H(zt(tt), xt(tt, :));
    end
    plot(t, Ht)
    xlabel('Time (s)')
    ylabel('H(x, z)')

    figure(fh + 5)
    hold on
    vt = zeros(size(t));
    for tt = 1:size(t, 1)
        vt(tt) = safety_filt_w_norm_dote(vd, xt(tt, :), zt(tt), sigma, alpha, delt);
    end
    plot(t, vt)
    ylim([-abs(vd), abs(vd)])
    drawnow
end