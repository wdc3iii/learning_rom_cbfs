classdef smoother_predictive_safety_filter
    %PREDICTIVE_SAFETY_FILTER Implements a predictive safety filter
    %   Properties:
    %   v_des: @(z) v_d(z) computes the desired vector field as a function
    %       of projection onto RoM
    %   k_track: @(x, v) tracking controller
    %   fx: @(x, u) FoM dynamics
    %   fz: @(z) RoM drift dynamics
    %   gz: @(z) RoM actuation matrix
    %   Pi: @(x) projects the FoM onto the RoM
    %   h: @(z) computes the value of the CBF (first order CBF on RoM)
    %   Jh: @(z) computes Jacobian of the CBF
    %   alpha: constant in dot_h >= - alpha h
    %   epsilon: barrier robustification term
    %   T: prediction horizon
    %   iters: maximum number of iterations to run
    %   tol: terminates algorithm when change in delta < tol
    %   K_delta: gain when adding violation to the CBF condition
    %   delta: state-dependent robustifcation term computed by algorithm

    properties
        % Controllers
        v_des
        k_track

        % Dynamics
        fx
        fz
        gz
        Pi
        JPi

        % Barrier conditions
        h
        Jh
        alpha
        epsilon

        % Algorithm parameters
        ts
        iters
        tol
        K_delta
        verbose
        fh

        % Running delta initial guess
        delta
        use_delta
    end

    methods
        function obj = smoother_predictive_safety_filter(v_des, k_track, fx, fz, gz, Pi, JPi, h, Jh, alpha, epsilon, T, dt, iters, tol, K_delta, verbose, fh, use_delta)
            %PREDICTIVE_SAFETY_FILTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.v_des = v_des;
            obj.k_track = k_track;
            obj.fx = fx;
            obj.fz = fz;
            obj.gz = gz;
            obj.Pi = Pi;
            obj.JPi = JPi;
            obj.h = h;
            obj.Jh = Jh;
            obj.alpha = alpha;
            obj.epsilon = epsilon;
            obj.iters = iters;
            obj.tol = tol;
            obj.K_delta = K_delta;
            obj.delta = 0;
            obj.verbose = verbose;
            obj.fh = fh;
            obj.ts = linspace(0, T, T / dt);
            obj.use_delta = use_delta;
        end

        function [obj, v] = apply(obj, t, x0)
            if obj.verbose > 0
                fprintf("t: %0.4f\n", t)
                fprintf("\tI: %d, delta: %0.4f\n", 0, obj.delta)
            end
            if obj.verbose > 1
                figure(obj.fh)
                clf
                hold on
            end
            if obj.use_delta
                for i = 1:obj.iters
                    % Roll out planner tracker with buffered safety function
                    [t, x] = ode45(@(t, x) obj.fx(x, obj.k_track(x, obj.safety_filter(obj.Pi(x)))), obj.ts, x0);

                    % Evaluate barrier function on full order model
                    h_bar = inf;
                    hx = zeros(size(t));
                    dot_hx = zeros(size(t));
                    for tt = 1:size(t, 1)
                        hx(tt) = obj.h(obj.Pi(x(tt, :)));
                        dot_hx(tt) = obj.Jh(obj.Pi(x(tt, :))) * obj.JPi(x(tt, :)) * obj.fx(x(tt, :), obj.k_track(x(tt, :), obj.safety_filter(obj.Pi(x(tt, :)))));
                        viol = dot_hx(tt) + obj.alpha * hx(tt);
                        if viol < h_bar
                            h_bar = viol;
                        end
                    end

                    % Buffer barrier function by violation of FoM
                    prev_zero = obj.delta == 0;
                    obj.delta = max(0, obj.delta - obj.K_delta * h_bar);
                    if obj.verbose > 0
                        fprintf("\tI: %d, delta: %0.4f\n", i, obj.delta)
                        % figure(2)
                        % clf
                        % hold on
                        % plot(dot_hx)
                        % plot(-obj.alpha * hx)
                        % pause
                    end
                    if obj.verbose > 1
                        l = (i - 1) / (obj.iters - 1);
                        c = l * [0.4940 0.1840 0.5560] + (1 - l) * [0.9290 0.6940 0.1250];
                        plot(t, hx, Color=c)
                    end
                    if abs(obj.K_delta * h_bar) < obj.tol || (prev_zero && obj.delta == 0)
                        break;
                    end
                end  % Hope for convergence
            end
            v = obj.safety_filter(obj.Pi(x0));
            if obj.verbose >= 10
                pause
            end
        end

        function v = safety_filter(obj, z)
            if obj.use_delta
                d = obj.delta;
            else
                d = 0;
            end
            vd = obj.v_des(z);
            Jh_ = obj.Jh(z);
            b = Jh_ * obj.gz(z);
            a = - obj.alpha * obj.h(z) - Jh_ * obj.fz(z) + d + obj.epsilon;
            if b * vd >= a
                v = vd;
            else
                v = vd + b' / (b * b') * (a - b * vd);
            end
        end

        function [obj, dot_x] = fx_cl(obj, t, x)
            [obj, v] = obj.apply(t, x);
            dot_x = obj.fx(x, obj.k_track(x, v));
        end

        function [t, x, h, vd, v, delta] = apply_cl(obj, x0, T, dt)
            t = linspace(0, T, T / dt + 1);
            x = zeros(size(t, 2), size(x0, 2));
            x(1, :) = x0;
            h = zeros(size(t));
            vd0 = obj.v_des(obj.Pi(x0));
            vd = zeros(size(t, 2), size(vd0, 1));
            v = zeros(size(t, 2), size(vd0, 1));
            delta = zeros(size(t));

            for tt = 1:size(t, 2) - 1
                h(tt) = obj.h(obj.Pi(x(tt, :)));
                [obj, dotx] = obj.fx_cl(t(tt), x(tt, :));
                delta(tt) = obj.delta;
                vd(tt, :) = obj.v_des(obj.Pi(x(tt, :)));
                v(tt, :) = obj.safety_filter(obj.Pi(x(tt, :)));
                x(tt + 1, :) = x(tt, :) + dotx' *  dt;
            end
            h(end) = obj.h(obj.Pi(x(end, :)));
            vd(end, :) = obj.v_des(obj.Pi(x(end, :)));
            v(end, :) = obj.safety_filter(obj.Pi(x(end, :)));
        end
    end
end

