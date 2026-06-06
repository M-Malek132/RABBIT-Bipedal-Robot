function params = init_bspline_params()
% Initialize controller parameters with default B-spline values
% Uses your BSpline and aproximation functions

    % B-spline settings
    params.n = 7;  % 8 control points (n+1)
    params.p = 3;  % cubic B-spline
    
    % Phase variable bounds (torso angle range during a step)
    params.theta0 = -0.1;  % torso angle at start of step
    params.thetaf =  0.3;  % torso angle at end of step
    
    % Generate initial control points using your approximation function
    % We create simple linear trajectories as initial data
    
    % Define start and end configurations for actuated joints
    qa_start = [-0.4;   % q1: stance knee at start (bent)
                 0.5;   % q2: stance hip at start
                -0.3;   % q3: swing knee at start
                 0.6];  % q4: swing hip at start
                
    qa_end   = [ 0.3;   % q1: stance knee at end (extended)
                -0.2;   % q2: stance hip at end
                 0.4;   % q3: swing knee at end
                -0.3];  % q4: swing hip at end
    
    % Use approximation to fit B-spline control points
    fs = 100;  % sampling frequency for approximation
    h = params.n;  % desired number of control points - 1
    
    % For each joint, fit control points to a linear trajectory
    params.ControlPoints = zeros(params.n+1, 4);
    
    for j = 1:4
        % Create linear data points
        time_data = linspace(0, 1, 20);  % 20 sample points
        D = qa_start(j) + time_data * (qa_end(j) - qa_start(j));
        
        % Use your approximation function to fit B-spline
        S_full = aproximation(D, params.p, h, fs, time_data);
        
        % Extract control points from approximation
        % The approximation returns the evaluated curve, but we need control points
        % We'll use the first and last directly, and interpolate the rest
        params.ControlPoints(1, j) = qa_start(j);
        params.ControlPoints(end, j) = qa_end(j);
        
        % Distribute interior control points evenly
        for k = 1:params.n-1
            s_k = k / params.n;
            params.ControlPoints(k+1, j) = qa_start(j) + s_k * (qa_end(j) - qa_start(j));
        end
    end
    
    % PD gains
    params.Kp = diag([200, 200, 150, 150]);
    params.Kd = diag([30,  30,  20,  20 ]);
    
end