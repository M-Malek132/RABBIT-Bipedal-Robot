classdef bezier_virtual_constraints < handle
    % Bézier polynomial virtual constraints for RABBIT bipedal robot
    % 
    % Phase variable: torso angle (qt = q(3))
    % s = (qt - theta0) / (thetaf - theta0)  ∈ [0,1]
    %
    % Virtual constraints:
    % y  = qa - hd(s)    where qa = actuated joints [q1;q2;q3;q4]
    % dy = dqa - dhd(s)
    
    properties
        M           % Order of Bézier polynomials
        n_outputs   % Number of virtual constraints (4 for RABBIT)
        theta0      % Torso angle at start of step
        thetaf      % Torso angle at end of step
    end
    
    methods
        function obj = bezier_virtual_constraints(M, theta0, thetaf)
            % Constructor
            % M: order of Bézier polynomials (5 => 6 coefficients)
            % theta0: initial torso angle at step start
            % thetaf: final torso angle at step end
            
            obj.M = M;
            obj.n_outputs = 4;  % q1, q2, q3, q4 (actuated joints)
            obj.theta0 = theta0;
            obj.thetaf = thetaf;
        end
        
        function [s, ds_dt] = compute_phase(obj, q, dq)
            % Compute normalized phase variable and its time derivative
            %
            % s = (qt - theta0) / (thetaf - theta0)
            % ds/dt = dqt / (thetaf - theta0)
            %
            % Inputs:
            %   q:  [x; z; qt; q1; q2; q3; q4]
            %   dq: [dx; dz; dqt; dq1; dq2; dq3; dq4]
            
            qt = q(3);
            dqt = dq(3);
            
            s = (qt - obj.theta0) / (obj.thetaf - obj.theta0);
            ds_dt = dqt / (obj.thetaf - obj.theta0);
        end
        
        function B = bernstein_basis(obj, s, derivative_order)
            % Evaluate Bernstein basis polynomials and their derivatives
            %
            % B(k+1) = M!/(k!(M-k)!) * s^k * (1-s)^(M-k)
            %
            % Inputs:
            %   s: normalized phase [0,1]
            %   derivative_order: 0 for basis, 1 for first derivative
            %
            % Output:
            %   B: 1 x (M+1) vector of basis values
            
            if nargin < 3
                derivative_order = 0;
            end
            
            M = obj.M;
            B = zeros(1, M+1);
            
            if derivative_order == 0
                % Zeroth order: Bernstein basis
                for k = 0:M
                    B(k+1) = nchoosek(M, k) * s^k * (1-s)^(M-k);
                end
                
            elseif derivative_order == 1
                % First derivative of Bernstein basis
                for k = 0:M
                    binom = nchoosek(M, k);
                    
                    if s == 0
                        % Handle s=0 specially
                        if k == 1
                            B(k+1) = binom * M;  % derivative at s=0
                        else
                            B(k+1) = 0;
                        end
                    elseif s == 1
                        % Handle s=1 specially
                        if k == M-1
                            B(k+1) = -binom * M;  % derivative at s=1
                        else
                            B(k+1) = 0;
                        end
                    else
                        % General case
                        term1 = 0;
                        term2 = 0;
                        if k > 0
                            term1 = k * s^(k-1) * (1-s)^(M-k);
                        end
                        if k < M
                            term2 = -(M-k) * s^k * (1-s)^(M-k-1);
                        end
                        B(k+1) = binom * (term1 + term2);
                    end
                end
            end
        end
        
        function hd = evaluate_desired_outputs(obj, alpha, s)
            % Evaluate desired joint angles at phase s
            %
            % hd(s) = sum_{k=0}^{M} B_k(s) * alpha(k+1)
            %
            % Inputs:
            %   alpha: (M+1) x 4 matrix of Bézier coefficients
            %          alpha(k+1, j) is coefficient k for output j
            %   s: phase variable value
            
            B = obj.bernstein_basis(s, 0);
            
            hd = zeros(obj.n_outputs, 1);
            for j = 1:obj.n_outputs
                hd(j) = B * alpha(:, j);
            end
        end
        
        function dhd = evaluate_desired_derivative(obj, alpha, s)
            % Evaluate derivative of desired outputs w.r.t. phase s
            %
            % dhd/ds = sum_{k=0}^{M} (dB_k/ds) * alpha(k+1)
            %
            % Note: This is dhd/ds, not dhd/dt
            % dhd/dt = dhd/ds * ds/dt
            
            dB = obj.bernstein_basis(s, 1);
            
            dhd = zeros(obj.n_outputs, 1);
            for j = 1:obj.n_outputs
                dhd(j) = dB * alpha(:, j);
            end
        end
        
        function [y, dy] = virtual_constraint(obj, q, dq, alpha)
            % Compute virtual constraints and their derivatives
            %
            % y  = qa - hd(s)
            % dy = dqa - dhd/ds * ds/dt
            %
            % Inputs:
            %   q:     [x; z; qt; q1; q2; q3; q4]
            %   dq:    [dx; dz; dqt; dq1; dq2; dq3; dq4]
            %   alpha: (M+1) x 4 Bézier coefficients
            
            % Compute phase variable
            [s, ds_dt] = obj.compute_phase(q, dq);
            
            % Get desired outputs and their phase derivative
            hd = obj.evaluate_desired_outputs(alpha, s);
            dhd_ds = obj.evaluate_desired_derivative(alpha, s);
            
            % Extract actuated joints
            qa  = q(4:7);
            dqa = dq(4:7);
            
            % Virtual constraints
            y  = qa  - hd;
            dy = dqa - dhd_ds * ds_dt;
        end
        
        function [y] = virtual_constraint_position(obj, q, alpha)
            % Compute just the position-level virtual constraint
            % Useful for optimization where velocities aren't needed
            
            [s, ~] = obj.compute_phase(q, zeros(7,1));
            hd = obj.evaluate_desired_outputs(alpha, s);
            qa = q(4:7);
            y = qa - hd;
        end
        
        function alpha = initialize_coefficients(obj, qa_start, qa_end)
            % Initialize Bézier coefficients for linear interpolation
            %
            % At s=0: hd = qa_start (B_0=1, all other B_k=0)
            % At s=1: hd = qa_end   (B_M=1, all other B_k=0)
            %
            % Inputs:
            %   qa_start: [q1;q2;q3;q4] at s=0 (beginning of step)
            %   qa_end:   [q1;q2;q3;q4] at s=1 (end of step, before impact)
            
            M = obj.M;
            alpha = zeros(M+1, obj.n_outputs);
            
            for j = 1:obj.n_outputs
                % Boundary coefficients
                alpha(1, j)   = qa_start(j);   % First control point
                alpha(M+1, j) = qa_end(j);     % Last control point
                
                % Interior coefficients: linear interpolation
                % This gives a reasonable but not optimal initial guess
                for k = 1:M-1
                    s_k = k / M;
                    alpha(k+1, j) = qa_start(j) + s_k * (qa_end(j) - qa_start(j));
                end
            end
        end
        
        function n = num_coefficients(obj)
            % Total number of optimization variables
            n = (obj.M + 1) * obj.n_outputs;
        end
        
        function alpha_vec = pack_coefficients(obj, alpha)
            % Convert matrix to vector for optimization
            alpha_vec = alpha(:);
        end
        
        function alpha = unpack_coefficients(obj, alpha_vec)
            % Convert vector back to matrix
            alpha = reshape(alpha_vec, obj.M+1, obj.n_outputs);
        end
        
        function plot_outputs(obj, alpha, num_points)
            % Plot desired trajectories vs phase variable
            if nargin < 3
                num_points = 100;
            end
            
            s_vec = linspace(0, 1, num_points);
            hd_plot = zeros(obj.n_outputs, num_points);
            
            for i = 1:num_points
                hd_plot(:, i) = obj.evaluate_desired_outputs(alpha, s_vec(i));
            end
            
            figure('Name', 'Virtual Constraints - Desired Trajectories');
            
            subplot(2,2,1);
            plot(s_vec, hd_plot(1,:), 'b-', 'LineWidth', 2);
            xlabel('Phase s'); ylabel('q_1 (rad)');
            title('Stance Knee (q1)'); grid on;
            
            subplot(2,2,2);
            plot(s_vec, hd_plot(2,:), 'r-', 'LineWidth', 2);
            xlabel('Phase s'); ylabel('q_2 (rad)');
            title('Stance Hip (q2)'); grid on;
            
            subplot(2,2,3);
            plot(s_vec, hd_plot(3,:), 'g-', 'LineWidth', 2);
            xlabel('Phase s'); ylabel('q_3 (rad)');
            title('Swing Knee (q3)'); grid on;
            
            subplot(2,2,4);
            plot(s_vec, hd_plot(4,:), 'm-', 'LineWidth', 2);
            xlabel('Phase s'); ylabel('q_4 (rad)');
            title('Swing Hip (q4)'); grid on;
            
            sgtitle('Bézier Virtual Constraints');
        end
        
        function plot_control_points(obj, alpha)
            % Visualize Bézier control points
            obj.plot_outputs(alpha);
            
            % Overlay control points
            M = obj.M;
            s_cp = linspace(0, 1, M+1);
            
            subplot(2,2,1); hold on;
            plot(s_cp, alpha(:,1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            
            subplot(2,2,2); hold on;
            plot(s_cp, alpha(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            
            subplot(2,2,3); hold on;
            plot(s_cp, alpha(:,3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            
            subplot(2,2,4); hold on;
            plot(s_cp, alpha(:,4), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
        end
    end
end