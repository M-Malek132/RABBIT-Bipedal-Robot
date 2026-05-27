function [u, delta_clf] = rabbit_clf_controller_v2(t, x, params)
% ---------------------------------------------------------------------
% RABBIT_CLF_CONTROLLER_V2
% یک کنترل‌کننده بهینه‌ساز QP بر پایه تابع لیاپانوف (CLF) برای ربات Rabbit
% ---------------------------------------------------------------------

%% ۱. استخراج حالت‌ها و پارامترها
p  = packParameters(params);
q  = x(1:7,1);
dq = x(8:14,1);

% محاسبه ماتریس‌های دینامیکی
D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_total = C + G; % مجموع نیروهای کوریولیس و گرانش
B = input_matrix();

%% ۲. خروجی‌ها و خطای ردیابی (Tracking Error)
% استخراج خروجی فعلی و مطلوب
[y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq);
[y_d, dy_d, ddy_d]  = get_DesiredOutputs(t);

% محاسبه بردار حالت خطا (eta)
e   = y - y_d;
de  = J_y * dq - dy_d;
eta = [e; de];

%% ۳. مشتقات لی (Lie Derivatives) و دینامیک خروجی
% سیستم به فرم: y_ddot = Lf_drift + Lg_ctrl * u
% برای پایداری عددی، به جای inv(D)*H از D \ H استفاده می‌کنیم
D_inv_H = D \ H_total;
D_inv_B = D \ B;

Lf_drift = Jdot_y_dq - J_y * D_inv_H - ddy_d;
Lg_ctrl  = J_y * D_inv_B;

%% ۴. ساخت تابع لیاپانوف (CLF)
% تعریف ماتریس‌های بهره برای تعیین رفتار همگرایی
n_y = length(e);
Kp_val = 100;
Kd_val = 20;
Kp_mat = Kp_val * eye(n_y);
Kd_mat = Kd_val * eye(n_y);

% حل تحلیلی معادله ریکاتی برای سیستم مرتبه دوم (P-matrix)
% V = eta' * P * eta
P = [Kp_mat + 0.5*(Kd_mat^2), 0.5*Kd_mat;
     0.5*Kd_mat,             0.5*eye(n_y)];

% سیستم خطا در فضای حالت: d_eta = F_drift + G_ctrl * u
F_drift = [de; Lf_drift];
G_ctrl  = [zeros(n_y, size(B,2)); Lg_ctrl];

% محاسبه مقدار لیاپانوف و مشتقات آن
V   = eta' * P * eta;
LfV = 2 * eta' * P * F_drift;
LgV = 2 * eta' * P * G_ctrl;

% نرخ همگرایی نمایی (Exponential Convergence Rate)
c_clf = 5.0; 

%% ۵. فرمولاسیون برنامه ریزی کوادراتیک (QP)
% متغیرهای تصمیم: z = [u (4x1); delta_clf (1x1)]
% هدف: Minimize 0.5 * z' * H_qp * z + f_qp' * z

n_u = size(B,2);
weight_u = 0.1;           % جریمه برای مصرف انرژی (گشتاور)
weight_delta_clf = 1e5;   % جریمه سنگین برای شل کردن قید پایداری (Relaxation)

H_qp = diag([weight_u * ones(n_u, 1); weight_delta_clf]);
f_qp = zeros(n_u + 1, 1);

% قید CLF: LfV + LgV*u <= -c_clf * V + delta_clf
% به فرم استاندارد: [LgV, -1] * z <= -LfV - c_clf * V
A_ineq = [LgV, -1];
b_ineq = -LfV - c_clf * V;

% محدودیت‌های اشباع گشتاور (Actuator Saturation)
u_max = 150; % حداکثر گشتاور مجاز موتورها [N-m]
lb = [-u_max * ones(n_u, 1); 0];    % delta_clf باید مثبت باشد
ub = [ u_max * ones(n_u, 1); Inf];

%% ۶. حل QP و خروجی نهایی
options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'active-set');
whos H_qp f_qp A_ineq b_ineq lb ub

[z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

if exitflag == 1
    u = z_opt(1:n_u);
    delta_clf = z_opt(n_u + 1);
else
    % استراتژی پشتیبان (Fallback) در صورت عدم حل QP
    warning('QP Solver failed at t = %.4f. Using PD Fallback.', t);
    u_pd = -(Kp_val*e + Kd_val*de);
    u = max(min(u_pd, u_max), -u_max);
    delta_clf = 0;
end

end

% ---------------------------------------------------------------------
% توابع کمکی داخلی
% ---------------------------------------------------------------------

function [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq)
    % در ربات Rabbit، معمولاً ۴ مفصل کنترل می‌شوند (به جز زاویه بدنه و موقعیت افقی/عمودی)
    % فرض بر این است که q = [x; y; theta; q1; q2; q3; q4]
    H_0 = [zeros(4, 3), eye(4)];
    y = H_0 * q;
    J_y = H_0;
    Jdot_y_dq = zeros(4, 1); % چون J_y ثابت است، مشتق آن صفر است
end

function [y_d, dy_d, ddy_d] = get_DesiredOutputs(t)
    % تنظیمات مسیر مطلوب برای هماهنگی با شرایط اولیه
    y_initial = [-0.3; 0.6; -1.0; 0.6]; 
    
    % پارامترهای نوسان
    Amp = [0.1; 0.15; 0.1; 0.15];
    freq = 1.5; % هرتز
    omega = 2 * pi * freq; 

    % تولید مسیر سینوسی
    y_d   = y_initial + Amp .* sin(omega * t);
    dy_d  = Amp .* omega .* cos(omega * t);
    ddy_d = -Amp .* (omega^2) .* sin(omega * t);
end
