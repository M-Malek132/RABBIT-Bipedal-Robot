function hzd_plotResults(t, x, u, model, opt)
%HZD_PLOTRESULTS  Visualise gait optimisation results.

params = model.params;
nq     = model.nq;
q      = x(:, 1:nq);
dq     = x(:, nq+1:end);

coordNames = {'p_x','p_z','q_t','q_1','q_2','q_3','q_4'};
uNames     = {'u_1','u_2','u_3','u_4'};

% ---- Joint angles ---
figure('Name','Joint Angles');
for i = 1:nq
    subplot(ceil(nq/2), 2, i);
    plot(t, q(:,i), 'LineWidth',1.5);
    grid on; xlabel('t [s]'); ylabel('[rad]');
    title(coordNames{i});
end
sgtitle('Generalised coordinates');

% ---- Torques ---
figure('Name','Torques');
for i = 1:model.nu
    subplot(2,2,i);
    plot(t, u(:,i), 'LineWidth',1.5);
    hold on;
    yline( opt.uMax,'r--','LineWidth',1);
    yline( opt.uMin,'r--','LineWidth',1);
    grid on; xlabel('t [s]'); ylabel('[Nm]');
    title(uNames{i});
end
sgtitle('Motor torques');

% ---- Foot and hip heights ---
N = length(t);
swingH = zeros(N,1);
hipH   = zeros(N,1);
for k = 1:N
    kin      = rabbit_kinematics(q(k,:)', params);
    swingH(k) = kin.swingFoot(2);
    hipH(k)   = kin.hip(2);
end

figure('Name','Heights');
plot(t, swingH, 'b', 'LineWidth',1.5); hold on;
plot(t, hipH,   'r', 'LineWidth',1.5);
yline(0, 'k--');
yline(opt.hipHeightMin,'r:','LineWidth',1);
grid on; xlabel('t [s]'); ylabel('height [m]');
legend('Swing foot','Hip','Ground','Min hip');
title('Swing foot & hip height');

% ---- Phase variable ---
figure('Name','Phase variable');
theta = zeros(N,1);
for k = 1:N
    theta(k) = hzd_phaseVariable(q(k,:)', model);
end
plot(t, theta, 'LineWidth',1.5);
yline(opt.thetaStart,'b--'); yline(opt.thetaEnd,'r--');
grid on; xlabel('t [s]'); ylabel('\theta');
legend('\theta','thetaStart','thetaEnd');
title('Phase variable (should be monotone per step)');
end
