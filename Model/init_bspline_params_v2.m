function ctrl = init_bspline_params_v2(x0, params)
% Initialize B-spline controller with better initial guess
% Uses simulation data or heuristics for better starting point

    % B-spline settings
    ctrl.n = 7;  % 8 control points
    ctrl.p = 3;  % cubic
    
    % Phase bounds based on initial state
    qt_start = x0(3);
    ctrl.theta0 = qt_start;       % start at current torso angle
    ctrl.thetaf = qt_start + 0.2; % end 0.2 rad forward (lean forward to walk)
    
    % Desired walking: robot should lean forward
    % Actuated joints should move to propel robot FORWARD
    
    % Stance leg: starts forward, pushes back
    % Swing leg: starts back, swings forward
    qa_start = [-0.3;   % q1: stance hip (slightly behind)
                 0.6;   % q2: stance knee (bent)
                -1.0;   % q3: swing hip (way back)
                 0.6];  % q4: swing knee (bent)
                
    qa_end   = [ 0.5;   % q1: stance hip (extended back - pushed off)
                -0.3;   % q2: stance knee (straightened)
                 0.3;   % q3: swing hip (forward - ready to land)
                -0.3];  % q4: swing knee (extended)
    
    % Initialize with linear interpolation
    ctrl.ControlPoints = zeros(ctrl.n+1, 4);
    for j = 1:4
        ctrl.ControlPoints(1, j) = qa_start(j);
        ctrl.ControlPoints(end, j) = qa_end(j);
        for k = 1:ctrl.n-1
            s_k = k / ctrl.n;
            ctrl.ControlPoints(k+1, j) = qa_start(j) + s_k * (qa_end(j) - qa_start(j));
        end
    end
    
    % Higher PD gains for better tracking
    ctrl.Kp = diag([400, 400, 300, 300]);
    ctrl.Kd = diag([50,  50,  35,  35 ]);
    
    fprintf('\nB-spline params v2:\n');
    fprintf('  Phase: [%.3f, %.3f]\n', ctrl.theta0, ctrl.thetaf);
    fprintf('  Start q_act: [%.2f, %.2f, %.2f, %.2f]\n', qa_start);
    fprintf('  End   q_act: [%.2f, %.2f, %.2f, %.2f]\n', qa_end);
    
end