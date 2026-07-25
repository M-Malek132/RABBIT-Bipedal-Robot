function animate_rabbit(x_traj, opts)
    % ANIMATE_RABBIT Visualizes the RABBIT trajectory, manages camera tracking,
    % renders environment/stepping stones, and generates a GIF of the simulation.
    %
    %   animate_rabbit(x_traj)          % 14xN or Nx14 trajectory
    %   animate_rabbit(x_traj, opts)
    %
    % opts is an optional struct; any subset of the fields below. The defaults
    % reproduce the historical single-argument behaviour exactly, so existing
    % callers (animate_hzd_result) are unaffected.
    %   .gif_file  output GIF path      (default Results/rabbit_animation.gif)
    %   .skip      draw every skip-th column of x_traj  (default 5)
    %   .delay     GIF frame delay [s]  (default 0.03)
    %   .title     axes title
    %   .name      figure window name

    if isempty(x_traj) || any(isnan(x_traj(:)))
        error('Animation aborted: trajectory is empty or contains NaNs.');
    end

    if nargin < 2 || isempty(opts), opts = struct(); end

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
    
    % Setup Result Folder. The default is anchored at the REPO ROOT, not the
    % cwd: 'Results' used to be a relative path tested with exist(...,'dir'),
    % which also searches the MATLAB path -- so from any cwd other than the repo
    % root that test passed on the repo's own Results/, mkdir was skipped, and
    % imwrite then failed with "Unable to open file ... for writing".
    % isfolder() looks only at the filesystem, so it cannot be fooled that way.
    default_gif = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                           'Results', 'rabbit_animation.gif');
    filename = getfielddef(opts, 'gif_file', default_gif);
    result_folder = fileparts(filename);
    if ~isempty(result_folder) && ~isfolder(result_folder)
        mkdir(result_folder);
    end

    % Figure Initialization
    figure('Color', 'w', 'Name', getfielddef(opts,'name','RABBIT Bipedal Walker'));
    clf;
    axis equal; grid on; hold on;
    xlabel('X [m]'); ylabel('Z [m]');
    title(getfielddef(opts, 'title', 'RABBIT 5-Link Walker Simulation'));
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
    
    skip_frames = max(1, round(getfielddef(opts, 'skip',  5)));
    delay_time  =        getfielddef(opts, 'delay', 0.03);
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
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        end
    end

    fprintf('Animation saved to:\n%s\n', filename);
end

function v = getfielddef(s, f, d)
    % Field of s if present and non-empty, else the default d.
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
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
    % Thin adapter around the SHARED kinematics helper so the body-point
    % geometry is defined in exactly one place (Utilities/get_body_points.m).
    % This used to duplicate the Tt/T2/T4/P_st/P_sw kinematics verbatim; the two
    % copies could silently drift, so this now just unpacks get_body_points.
    pts = get_body_points(q);
    stance_foot = pts.stance_foot;
    swing_foot  = pts.swing_foot;
    hip         = pts.hip;
    stance_knee = pts.stance_knee;
    swing_knee  = pts.swing_knee;
    torso_top   = pts.torso_top;
end
