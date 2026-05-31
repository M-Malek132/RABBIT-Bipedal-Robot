function gait_library = generate_gait_library()
% Generate gait library by recording working controller at different speeds
% Saves results to Results/ folder

    if exist('startup.m','file'), startup; end
    
    % Ensure Results folder exists
    results_dir = fullfile(pwd, 'Results');
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    params = parameters();
    
    % Speeds to generate gaits for
    speeds = [0.2, 0.3, 0.4, 0.5, 0.6];
    dqt_values = [0.15, 0.22, 0.30, 0.38, 0.45];
    
    gait_library = struct();
    valid_count = 0;
    
    for idx = 1:length(speeds)
        target_speed = speeds(idx);
        dqt0 = dqt_values(idx);
        
        fprintf('\n========================================\n');
        fprintf(' Generating gait for speed %.2f m/s\n', target_speed);
        fprintf('========================================\n');
        
        gait = record_one_gait(params, dqt0, target_speed);
        
        if ~isempty(gait) && isfield(gait, 'CP')
            valid_count = valid_count + 1;
            gait_library(valid_count).speed = gait.speed;
            gait_library(valid_count).CP = gait.CP;
            gait_library(valid_count).n = gait.n;
            gait_library(valid_count).p = gait.p;
            gait_library(valid_count).theta0 = gait.theta0;
            gait_library(valid_count).thetaf = gait.thetaf;
            gait_library(valid_count).T_step = gait.T_step;
            gait_library(valid_count).step_length = gait.step_length;
            gait_library(valid_count).actual_speed = gait.actual_speed;
            
            fprintf('  Step time: %.3f s\n', gait.T_step);
            fprintf('  Step length: %.3f m\n', gait.step_length);
            fprintf('  Actual speed: %.3f m/s\n', gait.actual_speed);
        else
            fprintf('  FAILED to record gait\n');
        end
    end
    
    if valid_count > 0
        % Save with timestamp
        timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
        save(fullfile(results_dir, ['gait_library_', timestamp, '.mat']), 'gait_library');
        
        % Also save as latest
        save(fullfile(results_dir, 'gait_library.mat'), 'gait_library');
        
        fprintf('\nGait library saved to Results/ with %d gaits.\n', valid_count);
    else
        fprintf('\nNo valid gaits generated.\n');
        gait_library = [];
    end
    
end

function gait = record_one_gait(params, dqt0, target_speed)
% Simulate one step and fit B-splines to the actuated joint trajectories

    nq = 7;
    x0 = zeros(2*nq, 1);
    
    qt = 0.1;
    q1 = -0.3;  q2 = 0.6;
    q3 = -1.0;  q4 = 0.6;
    
    l1 = params.l1;  l2 = params.l2;
    px = l1*sin(qt+q1) + l2*sin(qt+q1+q2);
    pz = l1*cos(qt+q1) + l2*cos(qt+q1+q2);
    
    q0 = [px; pz; qt; q1; q2; q3; q4];
    x0(1:nq) = q0;
    
    dq0 = zeros(nq, 1);
    dq0(3) = dqt0;
    
    J = J_stance(q0, packParameters(params));
    dq0_corrected = (eye(7) - pinv(J)*J) * dq0;
    x0(nq+1:end) = dq0_corrected;
    
    controller = @rabbit_controller;
    
    try
        [t_step, x_step, impact_info] = simulate_one_step(x0, params, controller);
        
        if isempty(t_step) || ~impact_info.detected
            fprintf('  No impact detected\n');
            gait = struct();
            return;
        end
        
        q_traj = x_step(:, 1:7);
        
        % Phase variable (same as working controller)
        theta_traj = q_traj(:, 4) + 0.5 * q_traj(:, 6);
        theta0 = theta_traj(1);
        thetaf = theta_traj(end);
        s_traj = (theta_traj - theta0) / (thetaf - theta0);
        
        % Actuated joint trajectories
        q_act = q_traj(:, 4:7);
        
        % B-spline settings
        n_cp = 8;
        p_deg = 3;
        n = n_cp - 1;
        
        % Fit control points by sampling at evenly spaced phase values
        CP = zeros(n_cp, 4);
        s_target = linspace(0, 1, n_cp)';
        
        for j = 1:4
            CP(:, j) = interp1(s_traj, q_act(:, j), s_target, 'linear', 'extrap');
        end
        
        % Verify fit
        s_check = linspace(0, 1, 50)';
        hd_check = zeros(50, 4);
        for k = 1:50
            N = BSpline(n, p_deg, s_check(k));
            hd_check(k, :) = (N * CP)';
        end
        
        q_act_interp = interp1(s_traj, q_act, s_check, 'linear', 'extrap');
        fit_error = rms(q_act_interp - hd_check, 'all');
        
        fprintf('  B-spline fit RMS error: %.4f rad\n', fit_error);
        
        gait.CP = CP;
        gait.n = n;
        gait.p = p_deg;
        gait.theta0 = theta0;
        gait.thetaf = thetaf;
        gait.T_step = t_step(end);
        gait.step_length = q_traj(end, 1) - q_traj(1, 1);
        gait.actual_speed = gait.step_length / gait.T_step;
        gait.speed = target_speed;
        gait.fit_error = fit_error;
        
    catch ME
        fprintf('  Error: %s\n', ME.message);
        gait = struct();
    end
    
end