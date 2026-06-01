% make_initial_guess_v2.m
% Sample desired_gait directly instead of recording

clear; clc;
if exist('startup.m','file'), startup; end

%% Use the SAME phase as the working controller
% rabbit_controller uses: theta = q1 + 0.5*q3
% With bounds: theta0 = -0.3, thetaf = 0.3

theta0 = -0.3;
thetaf = 0.3;

%% Sample desired_gait at evenly spaced phase values
n_cp = 8;
p_deg = 3;
n = n_cp - 1;

s_target = linspace(0, 1, n_cp)';
CP = zeros(n_cp, 4);

fprintf('Sampling desired_gait at %d phase values...\n', n_cp);
for k = 1:n_cp
    hd = desired_gait(s_target(k));
    CP(k, :) = hd';
    fprintf('  s=%.3f: hd=[%.3f, %.3f, %.3f, %.3f]\n', s_target(k), hd);
end

%% Save
save(fullfile('Results', 'initial_guess.mat'), ...
    'CP', 'theta0', 'thetaf', 'n', 'p_deg');
fprintf('\nSaved. Phase range: [%.1f, %.1f]\n', theta0, thetaf);

%% Verify by evaluating the B-spline and comparing to desired_gait
s_check = linspace(0, 1, 100)';
hd_original = zeros(100, 4);
hd_bspline = zeros(100, 4);

for k = 1:100
    hd_original(k, :) = desired_gait(s_check(k))';
    N = BSpline(n, p_deg, s_check(k));
    hd_bspline(k, :) = (N * CP)';
end

fit_err = rms(hd_original - hd_bspline, 'all');
fprintf('B-spline vs desired_gait RMS error: %.6f rad\n', fit_err);

%% Plot comparison
figure('Name', 'B-Spline vs desired_gait', 'Position', [50, 50, 1200, 400]);
joint_names = {'q1 (St Hip)', 'q2 (St Knee)', 'q3 (Sw Hip)', 'q4 (Sw Knee)'};

for j = 1:4
    subplot(1,4,j); hold on;
    plot(s_check, hd_original(:,j), 'b-', 'LineWidth', 2);
    plot(s_check, hd_bspline(:,j), 'r--', 'LineWidth', 2);
    plot(s_target, CP(:,j), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    xlabel('s'); ylabel('rad');
    title(joint_names{j}); grid on;
    if j == 1, legend('desired\_gait', 'B-spline', 'CPs'); end
end
sgtitle(sprintf('B-Spline Fit of desired\\_gait (err=%.6f rad)', fit_err));