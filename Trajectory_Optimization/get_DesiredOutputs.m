function [y_d, dy_d, ddy_d] = get_DesiredOutputs(t, x_initial_step, p)
    % x_initial_step is the state of the robot at the BEGINNING of the current step
    % t is the time elapsed SINCE the beginning of the step
    
    % 1. Get the robot's actual joint angles at the start of the step
    q0 = x_initial_step(1:7);
    y_start = q0(4:7); % These are the 4 actuated joints [q1; q2; q3; q4]
    
    % 2. Define where you WANT the joints to end up (target pose)
    % Example: a slight step forward
    y_final = [0.4; -0.2; -0.6; 0.2]; 
    
    % 3. Define the step duration
    T_step = 0.8; % seconds
    
    % 4. Use a smooth interpolation (C1 or C2 continuous)
    % A simple one is the 1-cos profile:
    if t < T_step
        % Smoothing factor s goes from 0 to 1
        % ds goes from 0 to 0 (at t=0 and t=T)
        s = 0.5 * (1 - cos(pi * t / T_step));
        ds = 0.5 * (pi / T_step) * sin(pi * t / T_step);
        dds = 0.5 * (pi / T_step)^2 * cos(pi * t / T_step);
    else
        s = 1; ds = 0; dds = 0;
    end
    
    % 5. Compute desired values
    % y_d starts exactly at y_start when t=0
    % dy_d starts exactly at 0 when t=0
    y_d   = y_start + (y_final - y_start) * s;
    dy_d  = (y_final - y_start) * ds;
    ddy_d = (y_final - y_start) * dds;
end
