function animate_rabbit_stepping_stones(x_traj)

if isempty(x_traj) || any(isnan(x_traj(:)))
    error('Animation aborted: trajectory is empty or contains NaNs.');
end

% Load configuration
cfg = config('get');

if isempty(cfg)
    error('Configuration not initialized. Run config(''init'') first.');
end

stones = cfg.stones;

x_all = x_traj';

figure('Color','w');
clf;
axis equal; grid on; hold on;
xlabel('X'); ylabel('Z');
title('RABBIT 5‑Link Walker with Stepping Stones');

%% Results folder
project_root = fileparts(mfilename('fullpath'));
result_folder = fullfile(project_root,'Results');

if ~exist(result_folder,'dir')
    mkdir(result_folder);
end

filename = fullfile(result_folder,'rabbit_animation.gif');

%% Draw stones
for i = 1:size(stones,1)
    xs = stones(i,1);
    xe = stones(i,2);

    patch([xs xe xe xs], ...
        [0 0 -0.1 -0.1], ...
        [0.6 0.6 0.6], ...
        'EdgeColor','k','LineWidth',1.5);
end

%% Graphics handles
h_stance_leg = plot(NaN,NaN,'b','LineWidth',4);
h_swing_leg  = plot(NaN,NaN,'r','LineWidth',4);
h_torso      = plot(NaN,NaN,'g','LineWidth',5);

h_hip         = plot(NaN,NaN,'ko','MarkerSize',10,'MarkerFaceColor','k');
h_stance_knee = plot(NaN,NaN,'ko','MarkerSize',8,'MarkerFaceColor','k');
h_swing_knee  = plot(NaN,NaN,'ko','MarkerSize',8,'MarkerFaceColor','k');
h_stance_foot = plot(NaN,NaN,'ks','MarkerSize',8,'MarkerFaceColor','b');
h_swing_foot  = plot(NaN,NaN,'ks','MarkerSize',8,'MarkerFaceColor','r');

%% Downsample
skip = 5;
frames = 1:skip:size(x_all,2);

%% Animation loop
for i = 1:length(frames)

    k = frames(i);
    q = x_all(1:7,k);
    px = q(1);

    update_rabbit_frame(q,handles);

    xlim([px-1 px+2]);
    ylim([-0.2 1.8]);

    drawnow limitrate;

    write_gif_frame(filename,i);

end

fprintf('Animation saved to:\n%s\n',filename);
end
