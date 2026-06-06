function ctrl = init_bspline_params()
% Initialize B-spline controller parameters
% Returns a structure with all controller-specific settings

    % B-spline settings
    ctrl.n = 7;  % 8 control points (n+1)
    ctrl.p = 3;  % cubic B-spline
    
    % Phase variable bounds (torso angle range during a step)
    % Based on your make_initial_state, qt starts around 0.1
    ctrl.theta0 = 0.05;  % torso angle at start of step
    ctrl.thetaf = 0.35;  % torso angle at end of step
    
    % Define start and end configurations for actuated joints
    % Based on your make_initial_state joint angles
    qa_start = [-0.3;   % q1: stance hip at start
                 0.6;   % q2: stance knee at start
                -1.0;   % q3: swing hip at start
                 0.6];  % q4: swing knee at start
                
    qa_end   = [ 0.3;   % q1: stance hip at end (extended back)
                -0.3;   % q2: stance knee at end (extended)
                 0.3;   % q3: swing hip at end (forward)
                -0.3];  % q4: swing knee at end (extended)
    
    % Initialize control points with linear interpolation
    ctrl.ControlPoints = zeros(ctrl.n+1, 4);
    
    for j = 1:4
        ctrl.ControlPoints(1, j) = qa_start(j);
        ctrl.ControlPoints(end, j) = qa_end(j);
        
        for k = 1:ctrl.n-1
            s_k = k / ctrl.n;
            ctrl.ControlPoints(k+1, j) = qa_start(j) + s_k * (qa_end(j) - qa_start(j));
        end
    end
    
    % PD gains (moderate values to start)
    ctrl.Kp = diag([300, 300, 200, 200]);
    ctrl.Kd = diag([40,  40,  25,  25 ]);
    
end