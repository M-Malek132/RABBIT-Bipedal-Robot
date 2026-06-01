classdef RabbitGaitOptimizer < handle
    % Optimize B-spline control points for periodic walking
    
    properties
        params
        controller
        target_speed
        x0
        options
        initial_cp   % Good initial guess
    end
    
    methods
        function obj = RabbitGaitOptimizer(params, controller, target_speed, initial_guess_file)
            obj.params = params;
            obj.controller = controller;
            obj.target_speed = target_speed;
            
            % Load good initial guess
            if nargin >= 4 && exist(initial_guess_file, 'file')
                data = load(initial_guess_file);
                obj.initial_cp = data.CP;
                obj.controller.trajectory.CP = data.CP;
                obj.controller.trajectory.theta0 = data.theta0;
                obj.controller.trajectory.thetaf = data.thetaf;
                fprintf('Loaded initial guess from %s (fit err: %.4f)\n', ...
                    initial_guess_file, data.fit_err);
            else
                % Use the working controller to generate initial guess
                obj.generate_initial_guess();
            end
            
            % Get initial state
            [obj.x0, ~, ~] = make_initial_state();
            
            % Optimization options - gentler for this problem
            obj.options = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...
                'MaxIterations', 50, ...
                'MaxFunctionEvaluations', 300, ...
                'OptimalityTolerance', 1e-2, ...
                'ConstraintTolerance', 5e-2, ...
                'StepTolerance', 1e-3, ...
                'FiniteDifferenceStepSize', 1e-3, ...
                'UseParallel', false);
        end
        
        function generate_initial_guess(obj)
            % Run working controller for one step and fit B-splines
            
            nq = 7;
            x0_local = zeros(2*nq, 1);
            
            qt = 0.1; q1 = -0.3; q2 = 0.6; q3 = -1.0; q4 = 0.6;
            l1 = obj.params.l1; l2 = obj.params.l2;
            px = l1*sin(qt+q1) + l2*sin(qt+q1+q2);
            pz = l1*cos(qt+q1) + l2*cos(qt+q1+q2);
            q0 = [px; pz; qt; q1; q2; q3; q4];
            x0_local(1:nq) = q0;
            
            dq0 = zeros(nq, 1);
            dq0(3) = 0.3;
            J = J_stance(q0, packParameters(obj.params));
            dq0_corrected = (eye(7) - pinv(J)*J) * dq0;
            x0_local(nq+1:end) = dq0_corrected;
            
            [t_step, x_step, ~] = simulate_one_step(x0_local, obj.params, @rabbit_controller);
            
            q_traj = x_step(:, 1:7);
            q_act = q_traj(:, 4:7);
            theta_traj = q_traj(:, 4) + 0.5 * q_traj(:, 6);
            theta0 = theta_traj(1);
            thetaf = theta_traj(end);
            s_traj = (theta_traj - theta0) / (thetaf - theta0);
            
            n_cp = 8;
            CP = zeros(n_cp, 4);
            s_target = linspace(0, 1, n_cp)';
            for j = 1:4
                CP(:, j) = interp1(s_traj, q_act(:, j), s_target, 'linear', 'extrap');
            end
            
            obj.initial_cp = CP;
            obj.controller.trajectory.CP = CP;
            obj.controller.trajectory.theta0 = theta0;
            obj.controller.trajectory.thetaf = thetaf;
            
            fprintf('Generated initial guess from working controller.\n');
            fprintf('  Phase range: [%.3f, %.3f]\n', theta0, thetaf);
        end
        
        function [CP_opt, fval, exitflag] = optimize(obj)
            traj = obj.controller.trajectory;
            
            % Start from good initial guess
            x0_opt = obj.initial_cp(:);
            
            % Tight bounds around initial guess (±0.2 rad per CP)
            lb = obj.initial_cp(:) - 0.2;
            ub = obj.initial_cp(:) + 0.2;
            
            % First and last CPs stay close to boundaries
            lb(1:4) = obj.initial_cp(1, :)' - 0.05;
            ub(1:4) = obj.initial_cp(1, :)' + 0.05;
            lb(end-3:end) = obj.initial_cp(end, :)' - 0.1;
            ub(end-3:end) = obj.initial_cp(end, :)' + 0.1;
            
            fprintf('\n============================================\n');
            fprintf('  OPTIMIZING PERIODIC GAIT\n');
            fprintf('  Target speed: %.2f m/s\n', obj.target_speed);
            fprintf('  Variables: %d\n', length(x0_opt));
            fprintf('============================================\n');
            
            [x_opt, fval, exitflag] = fmincon(...
                @(x) obj.cost_function(x), ...
                x0_opt, [], [], [], [], lb, ub, ...
                @(x) obj.nonlinear_constraints(x), ...
                obj.options);
            
            CP_opt = reshape(x_opt, traj.n+1, 4);
            traj.CP = CP_opt;
            
            fprintf('\nOptimization complete. Exit: %d, Cost: %.4f\n', exitflag, fval);
        end
    end
    
    methods (Access = private)
        
        function f = cost_function(obj, x)
            obj.controller.trajectory.set_from_optimization_vector(x);
            ctrl_handle = obj.controller.to_function_handle();
            
            try
                [t_step, x_step, impact_info] = simulate_one_step(...
                    obj.x0, obj.params, ctrl_handle);
                
                if isempty(t_step) || ~impact_info.detected || t_step(end) < 0.1
                    f = 1e6;
                    return;
                end
                
                % Simple cost: be close to initial guess + track speed
                cp_diff = sum((x - obj.initial_cp(:)).^2);
                
                dist = x_step(end,1) - obj.x0(1);
                speed = dist / t_step(end);
                speed_err = (speed - obj.target_speed)^2;
                
                f = 0.01 * cp_diff + speed_err;
                
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
                
                if isempty(t_step) || ~impact_info.detected || t_step(end) < 0.1
                    c = 100; ceq = 100 * ones(4, 1);
                    return;
                end
                
                x_plus = rabbit_reset_map(x_step(end,:)', obj.params);
                
                % Periodicity: just the 4 actuated joints (relaxed)
                ceq = [
                    x_plus(4) - obj.x0(6);    % q1+ = q3_0
                    x_plus(5) - obj.x0(7);    % q2+ = q4_0
                    x_plus(6) - obj.x0(4);    % q3+ = q1_0
                    x_plus(7) - obj.x0(5);    % q4+ = q2_0
                ];
                
                % Speed constraint
                dist = x_step(end,1) - obj.x0(1);
                speed = dist / t_step(end);
                c = abs(speed - obj.target_speed) - 0.2;
                
            catch
                c = 100; ceq = 100 * ones(4, 1);
            end
        end
    end
end