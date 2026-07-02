function animate_rabbit(x_traj)
    % ANIMATE_RABBIT Visualizes the RABBIT trajectory, manages camera tracking,
    % renders environment/stepping stones, and generates a GIF of the simulation.
    
    if isempty(x_traj) || any(isnan(x_traj(:)))
        error('Animation aborted: trajectory is empty or contains NaNs.');
    end

    %======================================================================
    % 1. Setup Environment & Configuration
    %======================================================================
    % Load configuration for stones (if configuration function exists)
    try
        cfg = config('get');
        stones = cfg.stones;
    catch
        stones = []; % Fallback if stepping stones config is not found
    end
    
    % Setup Result Folder
    result_folder = 'Results';
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
    end
    filename = fullfile(result_folder, 'rabbit_animation.gif');

    % Figure Initialization
    figure('Color', 'w', 'Name', 'RABBIT Bipedal Walker');
    clf;
    axis equal; grid on; hold on;
    xlabel('X [m]'); ylabel('Z [m]');
    title('RABBIT 5-Link Walker Simulation');
    view(2);

    %======================================================================
    % 2. Initialize Environment Geometry
    %======================================================================
    % Draw Ground or Stepping Stones
    if isempty(stones)
        plot([-100 100], [0 0], 'k', 'LineWidth', 2);
    else
        for i = 1:size(stones, 1)
            patch([stones(i,1) stones(i,2) stones(i,2) stones(i,1)], ...
                  [0 0 -0.1 -0.1], [0.6 0.6 0.6], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end
    end

    %======================================================================
    % 3. Initialize Graphics Handles
    %======================================================================
    % Robot Links
    handles.stance_leg = plot(NaN, NaN, 'b', 'LineWidth', 4);
    handles.swing_leg  = plot(NaN, NaN, 'r', 'LineWidth', 4);
    handles.torso      = plot(NaN, NaN, 'g', 'LineWidth', 5);
    
    % Robot Joints
    handles.hip         = plot(NaN, NaN, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    handles.stance_knee = plot(NaN, NaN, 'ko', 'MarkerSize', 8,  'MarkerFaceColor', 'k');
    handles.swing_knee  = plot(NaN, NaN, 'ko', 'MarkerSize', 8,  'MarkerFaceColor', 'k');
    handles.stance_foot = plot(NaN, NaN, 'ks', 'MarkerSize', 8,  'MarkerFaceColor', 'b');
    handles.swing_foot  = plot(NaN, NaN, 'ks', 'MarkerSize', 8,  'MarkerFaceColor', 'r');

    %======================================================================
    % 4. Animation and GIF Generation Loop
    %======================================================================
    % Transpose trajectory if it is oriented vertically (states as rows)
    if size(x_traj, 1) ~= 14
        x_traj = x_traj';
    end
    
    skip_frames = 5; 
    frame_indices = 1 : skip_frames : size(x_traj, 2);
    
    for idx = 1:length(frame_indices)
        k = frame_indices(idx);
        q = x_traj(1:7, k);
        px = q(1);
        
        % Update Link & Joint Graphics positions
        update_rabbit_frame(q, handles);
        
        % Dynamic Camera Tracking
        xlim([px - 1.5, px + 1.5]);
        ylim([-0.2, 1.8]);
        drawnow limitrate;
        
        % Capture and save GIF frame
        frame = getframe(gcf);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);
        
        if idx == 1
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.03);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
        end
    end
    
    fprintf('Animation saved to:\n%s\n', filename);
end

%==========================================================================
% Helper Functions
%==========================================================================
function update_rabbit_frame(q, handles)
    % Retrieves points and updates plotting coordinates for the robot body
    [stance_foot, swing_foot, hip, stance_knee, swing_knee, torso_top] = rabbit_points(q);
    
    % Update Stance Leg (hip -> knee -> foot)
    set(handles.stance_leg, 'XData', [hip(1), stance_knee(1), stance_foot(1)], ...
                            'YData', [hip(2), stance_knee(2), stance_foot(2)]);
                        
    % Update Swing Leg (hip -> knee -> foot)
    set(handles.swing_leg,  'XData', [hip(1), swing_knee(1), swing_foot(1)], ...
                            'YData', [hip(2), swing_knee(2), swing_foot(2)]);
                        
    % Update Torso
    set(handles.torso,      'XData', [hip(1), torso_top(1)], ...
                            'YData', [hip(2), torso_top(2)]);
                        
    % Update Joints & End-Effectors
    set(handles.hip,         'XData', hip(1),         'YData', hip(2));
    set(handles.stance_knee, 'XData', stance_knee(1), 'YData', stance_knee(2));
    set(handles.swing_knee,  'XData', swing_knee(1),  'YData', swing_knee(2));
    set(handles.stance_foot, 'XData', stance_foot(1), 'YData', stance_foot(2));
    set(handles.swing_foot,  'XData', swing_foot(1),  'YData', swing_foot(2));
end

function [stance_foot, swing_foot, hip, stance_knee, swing_knee, torso_top] = rabbit_points(q)
    % Wrapper for kinematic position functions generated from your dynamic model
    stance_foot = P_st(q);
    swing_foot  = P_sw(q);

    pos_hip = Tt(q)* [0 0 0 1]';
    P_torso = Tt(q)* [0 -0.75 0 1]';
    P_knee_st = T2(q)* [0 0 0 1]';
    P_knee_sw = T4(q)* [0 0 0 1]';

    hip         = [pos_hip(1,1);   pos_hip(3,1)];
    stance_knee = [P_knee_st(1,1);   P_knee_st(3,1)];
    swing_knee  = [P_knee_sw(1,1);   P_knee_sw(3,1)];
    torso_top   = [P_torso(1,1);   P_torso(3,1)];
end
