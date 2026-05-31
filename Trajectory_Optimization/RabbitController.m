classdef RabbitController < handle
    % PD controller tracking B-spline virtual constraints for RABBIT
    %
    % Uses BSplineTrajectory to define desired joint angles as
    % functions of the phase variable (torso angle).
    
    properties
        trajectory    % BSplineTrajectory object
        Kp            % 4x4 proportional gain matrix
        Kd            % 4x4 derivative gain matrix
    end
    
    methods
        function obj = RabbitController(trajectory, Kp, Kd)
            % Constructor
            obj.trajectory = trajectory;
            obj.Kp = Kp;
            obj.Kd = Kd;
        end
        
        function u = compute(obj, t, x)
            % Compute control torques
            % t: current time
            % x: full state [q; dq]
            % Returns 4x1 torque vector
            
            nq = 7;
            q  = x(1:nq);
            dq = x(nq+1:end);
            
            % Compute virtual constraints
            [y, dy] = obj.trajectory.virtual_constraint(q, dq);
            
            % PD control
            u = -obj.Kp * y - obj.Kd * dy;
        end
        
        function handle = to_function_handle(obj)
            % Convert to function handle for simulator
            % simulator expects: u = controller(t, x, params)
            handle = @(t, x, ~) obj.compute(t, x);
        end
    end
    
    methods (Static)
        function ctrl = default()
            % Create controller with default parameters
            traj = BSplineTrajectory(7, 3, 0.05, 0.30);
            
            % Default start/end configurations
            q_start = [-0.3; 0.6; -1.0; 0.6];
            q_end   = [ 0.5; -0.3; 0.3; -0.3];
            traj.initialize_linear(q_start, q_end);
            
            Kp = diag([400, 400, 300, 300]);
            Kd = diag([50,  50,  35,  35]);
            
            ctrl = RabbitController(traj, Kp, Kd);
        end
    end
end