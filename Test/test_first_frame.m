clc;clear 
close all;

[x0, params, p] = make_initial_state();

debug_first_frame(x0, params);

function debug_first_frame(x0, params)
    figure('Color','w');
    axis equal; grid on; hold on;
    xlabel('X'); ylabel('Z');
    title('RABBIT First Frame Debug');
    view(2);

    if isfield(params, 'stones')
        for i = 1:size(params.stones, 1)
            x_s = params.stones(i, 1);
            x_e = params.stones(i, 2);
            patch([x_s, x_e, x_e, x_s], [0, 0, -0.1, -0.1], [0.6 0.6 0.6], ...
                'EdgeColor', 'k', 'LineWidth', 1.5);
        end
    else
        plot([-100 100], [0 0], 'k', 'LineWidth', 2);
    end

    p = packParameters(params);

    q = x0(1:7);

    [stance_foot, swing_foot, hip, stance_knee, swing_knee, torso_top] = rabbit_kinematics(q, p);

    % Draw robot
    plot([hip(1), stance_knee(1), stance_foot(1)], ...
         [hip(2), stance_knee(2), stance_foot(2)], 'b', 'LineWidth', 4);

    plot([hip(1), swing_knee(1), swing_foot(1)], ...
         [hip(2), swing_knee(2), swing_foot(2)], 'r', 'LineWidth', 4);

    plot([hip(1), torso_top(1)], [hip(2), torso_top(2)], 'g', 'LineWidth', 5);

    plot(hip(1), hip(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    plot(stance_knee(1), stance_knee(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(swing_knee(1), swing_knee(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(stance_foot(1), stance_foot(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot(swing_foot(1), swing_foot(2), 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

    xlim([x0(1) - 1.0, x0(1) + 2.0]);
    ylim([-0.2, 1.8]);

    fprintf('x0 = \n');
    disp(x0);
    fprintf('stance foot = [%f, %f]\n', stance_foot(1), stance_foot(2));
    fprintf('swing  foot = [%f, %f]\n', swing_foot(1), swing_foot(2));
    fprintf('hip         = [%f, %f]\n', hip(1), hip(2));
end
