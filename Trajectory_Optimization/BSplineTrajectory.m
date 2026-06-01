classdef BSplineTrajectory < handle
    % B-spline virtual constraints for RABBIT bipedal robot
    %
    % Phase variable: theta = q1 + 0.5*q3  (stance hip + 0.5*swing hip)
    % Normalized: s = (theta - theta0) / (thetaf - theta0)
    %
    % Virtual constraint: y = q_actuated - h_d(s)
    
    properties
        n           % number of data points (CPs = n+1)
        p           % B-spline degree
        theta0      % phase variable at start of step
        thetaf      % phase variable at end of step
        CP          % (n+1) x 4 control points matrix
    end
    
    methods
        function obj = BSplineTrajectory(n, p, theta0, thetaf)
            obj.n = n;
            obj.p = p;
            obj.theta0 = theta0;
            obj.thetaf = thetaf;
            obj.CP = zeros(n+1, 4);
        end
        
        function theta = phase_variable(obj, q)
            % Compute raw phase variable (same as working controller)
            % theta = q1 + 0.5*q3 = q(4) + 0.5*q(6)
            theta = q(4) + 0.5 * q(6);
        end
        
        function s = phase(obj, q)
            % Compute normalized phase variable s ∈ [0,1]
            theta = obj.phase_variable(q);
            s = (theta - obj.theta0) / (obj.thetaf - obj.theta0);
            s = min(max(s, 0), 1);
        end
        
        function ds = phase_derivative(obj, q, dq)
            % Compute time derivative of phase variable
            % dtheta/dt = dq1 + 0.5*dq3 = dq(4) + 0.5*dq(6)
            dtheta = dq(4) + 0.5 * dq(6);
            ds = dtheta / (obj.thetaf - obj.theta0);
        end
        
        function hd = evaluate(obj, s)
            % Evaluate desired joint angles at phase s
            N = BSpline(obj.n, obj.p, s);
            hd = obj.CP' * N';
        end
        
        function dhd = evaluate_derivative(obj, s)
            % Evaluate derivative of desired angles w.r.t. phase s
            dN = BSpline_derivative(obj.n, obj.p, s);
            dhd = obj.CP' * dN';
        end
        
        function [y, dy] = virtual_constraint(obj, q, dq)
            % Compute virtual constraint and its derivative
            s = obj.phase(q);
            ds = obj.phase_derivative(q, dq);
            
            hd = obj.evaluate(s);
            dhd_ds = obj.evaluate_derivative(s);
            
            q_act = q(4:7);
            dq_act = dq(4:7);
            
            y  = q_act - hd;
            dy = dq_act - dhd_ds * ds;
        end
        
        function initialize_linear(obj, q_start, q_end)
            % Initialize control points with linear interpolation
            for j = 1:4
                obj.CP(1, j) = q_start(j);
                obj.CP(obj.n+1, j) = q_end(j);
                for k = 1:obj.n-1
                    s_k = k / obj.n;
                    obj.CP(k+1, j) = q_start(j) + s_k*(q_end(j) - q_start(j));
                end
            end
        end
        
        function x = get_optimization_vector(obj)
            x = obj.CP(:);
        end
        
        function set_from_optimization_vector(obj, x)
            obj.CP = reshape(x, obj.n+1, 4);
        end
        
        function n_vars = num_variables(obj)
            n_vars = (obj.n + 1) * 4;
        end
        
        function plot(obj)
            s_vec = linspace(0, 1, 200);
            hd = zeros(4, 200);
            for i = 1:200
                hd(:, i) = obj.evaluate(s_vec(i));
            end
            
            joint_names = {'Stance Hip (q1)', 'Stance Knee (q2)', ...
                           'Swing Hip (q3)', 'Swing Knee (q4)'};
            
            figure('Name', 'B-Spline Virtual Constraints');
            for j = 1:4
                subplot(2, 2, j);
                plot(s_vec, hd(j, :), 'b-', 'LineWidth', 2); hold on;
                cp_s = linspace(0, 1, obj.n+1);
                plot(cp_s, obj.CP(:, j), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                xlabel('Phase s'); ylabel('Angle (rad)');
                title(joint_names{j}); grid on;
            end
        end
    end
end