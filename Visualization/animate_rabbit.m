function animate_rabbit(x_traj, params)
    
    figure('Color','w');
    axis equal; grid on; hold on;
    xlabel('X'); ylabel('Z');
    title('RABBIT 5-Link Walker');
    view(2);

<<<<<<< HEAD
    filename = 'Results/rabbit_animation.gif';
=======
    % Create 'result' folder if it doesn't exist
    result_folder = 'Results';
    if ~exist(result_folder, 'dir')
        mkdir(result_folder);
        fprintf('Created folder: %s\n', result_folder);
    end

    % Save GIF in the result folder
    filename = fullfile(result_folder, 'rabbit_animation.gif');

>>>>>>> Trajectory-Optimization

    %==============================
    % Initialize Plot Handles
    %==============================
    % Ground
    plot([-100 100], [0 0], 'k', 'LineWidth', 2);
    
    % Initialize empty lines for robot parts
    h_stance_leg = plot(NaN, NaN, 'b', 'LineWidth', 4);
    h_swing_leg  = plot(NaN, NaN, 'r', 'LineWidth', 4);
    h_torso      = plot(NaN, NaN, 'g', 'LineWidth', 5);
    
    % Initialize joints
    h_hip         = plot(NaN, NaN, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    h_stance_knee = plot(NaN, NaN, 'ko', 'MarkerSize', 8,  'MarkerFaceColor', 'k');
    h_swing_knee  = plot(NaN, NaN, 'ko', 'MarkerSize', 8,  'MarkerFaceColor', 'k');
    h_stance_foot = plot(NaN, NaN, 'ks', 'MarkerSize', 8,  'MarkerFaceColor', 'b');
    h_swing_foot  = plot(NaN, NaN, 'ks', 'MarkerSize', 8,  'MarkerFaceColor', 'r');

    p = packParameters(params);
    
    % Downsample data to speed up animation and GIF generation
    % Adjust this number if it's still too slow or too fast
    skip_frames = 10; 
    frame_indices = 1 : skip_frames : size(x_traj, 2);
    
    for idx = 1:length(frame_indices)
        k = frame_indices(idx);
        
        % States & Kinematics
        q = x_traj(1:7, k);
        px = q(1);
        [stance_foot, swing_foot, hip, stance_knee, swing_knee, torso_top] = rabbit_kinematics(q, p);

        %==============================
        % Update Plot Data Efficiently
        %==============================
        set(h_stance_leg, 'XData', [hip(1), stance_knee(1), stance_foot(1)], ...
                          'YData', [hip(2), stance_knee(2), stance_foot(2)]);
                      
        set(h_swing_leg,  'XData', [hip(1), swing_knee(1), swing_foot(1)], ...
                          'YData', [hip(2), swing_knee(2), swing_foot(2)]);
                      
        set(h_torso,      'XData', [hip(1), torso_top(1)], ...
                          'YData', [hip(2), torso_top(2)]);

        set(h_hip,         'XData', hip(1),         'YData', hip(2));
        set(h_stance_knee, 'XData', stance_knee(1), 'YData', stance_knee(2));
        set(h_swing_knee,  'XData', swing_knee(1),  'YData', swing_knee(2));
        set(h_stance_foot, 'XData', stance_foot(1), 'YData', stance_foot(2));
        set(h_swing_foot,  'XData', swing_foot(1),  'YData', swing_foot(2));

        % Camera tracking
        xlim([px - 1.5, px + 1.5]);
        ylim([-0.1, 1.8]);

        drawnow limitrate;

        % Capture frame for GIF
        frame = getframe(gcf);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);

        % Write GIF
        if idx == 1
            imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.03);
        else
            imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.03);
        end
    end
    
    fprintf('Animation saved to: %s\n', filename);
end
