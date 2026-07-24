function ctrl = init_bspline_params()
    % --- Spline Structure ---
    ctrl.n_segments = 7; % Number of segments
    ctrl.degree = 3;      % Cubic B-spline
    ctrl.n_ctrl_points = ctrl.n_segments + ctrl.degree; 
    ctrl.nu = 4;           % Actuated joints
    
    % --- Phase Variable Range (e.g., Hip Angle) ---
    % This defines the "clock" of the gait
    ctrl.s_min = 0.05; 
    ctrl.s_max = 0.35;
    
    % --- Boundary Conditions (Initial/Final Joint States) ---
    % [Joint 1; Joint 2; Joint 3; Joint 4]
    qa_start = [-0.3;  0.6; -1.0;  0.6];
    qa_end   = [ 0.3; -0.3;  0.3; -0.3];
    
    % --- Initialize Control Points ---
    % Shape: (Number of Control Points) x (Number of Joints)
    ctrl.ControlPoints = zeros(ctrl.n_ctrl_points, ctrl.nu);
    
    % Fill control points
    for j = 1:ctrl.nu
        % Set the first and last control points to satisfy boundary conditions
        % Note: For B-splines, the actual curve doesn't hit the first/last 
        % control point unless you use specific clamping. 
        % For an initial guess, a linear interpolation is a safe way to start.
        
        lin_space = linspace(qa_start(j), qa_end(j), ctrl.n_ctrl_points);
        ctrl.ControlPoints(:, j) = lin_space';
    end
    
    % --- PD Controller Gains ---
    % Using diagonal matrices for easy computation in the ODE
    ctrl.Kp = diag([300, 300, 200, 200]);
    ctrl.Kd = diag([40,  40,  25,  25]);
end
