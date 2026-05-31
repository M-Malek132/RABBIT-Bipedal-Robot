classdef RabbitGaitOptimizer < handle
    % Optimize B-spline control points for periodic walking
    
    properties
        params
        controller
        target_speed
        x0          % Cached initial state
        options
    end
    
    methods
        function obj = RabbitGaitOptimizer(params, controller, target_speed)
            obj.params = params;
            obj.controller = controller;
            obj.target_speed = target_speed;
            
            % Cache initial state ONCE
            [obj.x0, ~, ~] = make_initial_state();
            
            % Adjust phase bounds once
            obj.controller.trajectory.theta0 = obj.x0(3);
            obj.controller.trajectory.thetaf = obj.x0(3) + 0.2;
            
            obj.options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...
                'MaxIterations', 100, ...
                'MaxFunctionEvaluations', 500, ...
                'OptimalityTolerance', 1e-3, ...
                'ConstraintTolerance', 1e-2, ...
                'StepTolerance', 1e-4, ...
                'UseParallel', false);
        end
        
        function [CP_opt, fval, exitflag] = optimize(obj)
            traj = obj.controller.trajectory;
            
            % Initial guess
            x0_opt = traj.get_optimization_vector();
            
            % Bounds
            lb = repmat([-pi/2; -pi/2; -pi; -pi/2], traj.n+1, 1);
            ub = repmat([ pi/2;  pi/2;  pi;  pi/2], traj.n+1, 1);
            
            q_start = traj.CP(1, :)';
            q_end = traj.CP(end, :)';
            lb(1:4) = q_start - 0.2;
            ub(1:4) = q_start + 0.2;
            lb(end-3:end) = q_end - 0.2;
            ub(end-3:end) = q_end + 0.2;
            
            % Run
            [x_opt, fval, exitflag] = fmincon(...
                @(x) obj.cost_function(x), ...
                x0_opt, [], [], [], [], lb(:), ub(:), ...
                @(x) obj.nonlinear_constraints(x), ...
                obj.options);
            
            CP_opt = reshape(x_opt, traj.n+1, 4);
            traj.CP = CP_opt;
            
            fprintf('\nOptimization complete. Exit flag: %d\n', exitflag);
        end
    end
    
    methods (Access = private)
        
        function f = cost_function(obj, x)
            % Simple cost: just minimize periodicity error
            % Speed is handled as a constraint, not in cost
            
            obj.controller.trajectory.set_from_optimization_vector(x);
            
            ctrl_handle = obj.controller.to_function_handle();
            
            try
                [t_step, ~, impact_info] = simulate_one_step(...
                    obj.x0, obj.params, ctrl_handle);
                
                if isempty(t_step) || ~impact_info.detected
                    f = 1e6;
                    return;
                end
                
                % Simple cost: prefer longer steps (more efficient)
                f = -t_step(end);
                
            catch
                f = 1e6;
            end
        end
        
        function [c, ceq] = nonlinear_constraints(obj, x)
            
            obj.controller.trajectory.set_from_optimization_vector(x);
            
            ctrl_handle = obj.controller.to_function_handle();
            
            try
                [t_step, x_step, impact_info] = simulate_one_step(...
                    obj.x0, obj.params, ctrl_handle);
                
                if isempty(t_step) || ~impact_info.detected
                    c = 100;
                    ceq = 100 * ones(14, 1);
                    return;
                end
                
                % Reset map
                x_plus = rabbit_reset_map(x_step(end,:)', obj.params);
                
                % Periodicity constraints (with leg swap)
                ceq = [
                    x_plus(4) - obj.x0(6);    % q1+ = q3_0
                    x_plus(5) - obj.x0(7);    % q2+ = q4_0
                    x_plus(6) - obj.x0(4);    % q3+ = q1_0
                    x_plus(7) - obj.x0(5);    % q4+ = q2_0
                    x_plus(11) - obj.x0(13);  % dq1+ = dq3_0
                    x_plus(12) - obj.x0(14);  % dq2+ = dq4_0
                    x_plus(13) - obj.x0(11);  % dq3+ = dq1_0
                    x_plus(14) - obj.x0(12);  % dq4+ = dq2_0
                ];
                
                % Speed constraint: match target speed
                dist = x_step(end,1) - obj.x0(1);
                speed = dist / t_step(end);
                c_speed = abs(speed - obj.target_speed) - 0.1;  % within 0.1 m/s
                
                % Foot clearance at mid-step
                mid = round(length(t_step)/2);
                q_mid = x_step(mid, 1:7)';
                [~, p_sw, ~, ~, ~, ~] = rabbit_kinematics(q_mid, packParameters(obj.params));
                c_foot = 0.01 - p_sw(2);
                
                c = [c_speed; c_foot];
                
            catch
                c = 100;
                ceq = 100 * ones(8, 1);
            end
        end
    end
end