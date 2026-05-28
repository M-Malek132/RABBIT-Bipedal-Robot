function [u, delta_clf] = rabbit_clf_controller_v2(t, x, params)
% ---------------------------------------------------------------------
% RABBIT_CLF_CONTROLLER_V3 (Robust Version)
% کنترل‌کننده CLF-QP برای ربات Rabbit با مکانیزم جلوگیری از Infeasibility
% ---------------------------------------------------------------------

%% ۱. استخراج حالت‌ها و پارامترها
p  = packParameters(params);
q  = x(1:7,1);
dq = x(8:14,1);

% محاسبه ماتریس‌های دینامیکی (استفاده از عملگر \ به جای inv)
D = D_matrix(q,p);
C = C_vector(q,dq,p);
G = G_vector(q,p);
H_total = C + G; 
B = input_matrix();

%% ۲. خروجی‌ها و خطای ردیابی
[y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq);
[y_d, dy_d, ddy_d]  = get_DesiredOutputs(t,x,p);

% بردار خطا و مشتق خطا
e   = y - y_d;
de  = J_y * dq - dy_d;
eta = [e; de];

%% ۳. مشتقات لی (Lie Derivatives)
% محاسبه شتاب‌ها با روش تقسیم ماتریسی (پایدارتر از معکوس ماتریس)
D_inv_H = D \ H_total;
D_inv_B = D \ B;

% Lf_drift: بخشی از شتاب که به ورودی u وابسته نیست
Lf_drift = Jdot_y_dq - J_y * D_inv_H - ddy_d;
% Lg_ctrl: ماتریس اثر ورودی بر شتاب خروجی
Lg_ctrl  = J_y * D_inv_B;

%% ۴. ساخت تابع لیاپانوف (CLF)
n_y = length(e);
% ماتریس‌های بهره برای ساخت P (قابل تنظیم برای سرعت پاسخ‌دهی)
Kp_val = 100; 
Kd_val = 20;
Kp_mat = Kp_val * eye(n_y);
Kd_mat = Kd_val * eye(n_y);

% ماتریس P از حل معادله لیاپانوف برای سیستم خطا
P = [Kp_mat + 0.5*(Kd_mat^2), 0.5*Kd_mat;
     0.5*Kd_mat,             0.5*eye(n_y)];

% سیستم خطا در فضای حالت: d_eta = F_drift + G_ctrl * u
F_drift = [de; Lf_drift];
G_ctrl  = [zeros(n_y, size(B,2)); Lg_ctrl];

% محاسبه لیاپانوف و مشتقات آن
V   = eta' * P * eta;
LfV = 2 * eta' * P * F_drift;
LgV = 2 * eta' * P * G_ctrl;

% نرخ همگرایی نمایی (Exponential Convergence Rate)
% اگر t کوچک است، نرخ را کمتر می‌گیریم تا سیستم به آرامی به مسیر جذب شود
c_clf = 2.0; 
if t < 0.2, c_clf = 0.5; end 

%% ۵. فرمولاسیون QP (رفع مشکل Infeasibility)
n_u = size(B,2);
n_z = n_u + 1; % u1, u2, u3, u4 + delta_clf

% وزن‌دهی در تابع هدف
% وزن u را بسیار کوچک می‌گیریم تا حل‌کننده آزادی عمل داشته باشد
% وزن delta_clf باید بزرگ باشد اما نه آنقدر که ماتریس بد-حالت شود
weight_u = 0.01;           
weight_delta_clf = 1000;   

H_qp = diag([weight_u * ones(1, n_u), weight_delta_clf]);
% برای اطمینان از مثبت معین بودن ماتریس H
H_qp = H_qp + eye(n_z) * 1e-8; 

f_qp = zeros(n_z, 1);

% قید CLF: LfV + LgV*u <= -c_clf * V + delta_clf
% فرم استاندارد: A*z <= b
A_ineq = [LgV, -1]; 
b_ineq = -LfV - c_clf * V;

% محدودیت‌های فیزیکی گشتاور (Actuator Limits)
u_max = 150; 
lb = [-u_max * ones(n_u, 1); 0];    % delta_clf همواره >= 0
ub = [ u_max * ones(n_u, 1); Inf];

%% ۶. حل QP
options = optimoptions('quadprog', ...
    'Algorithm', 'interior-point-convex', ...
    'Display', 'off', ...
    'ConstraintTolerance', 1e-4, ...
    'OptimalityTolerance', 1e-4);

[z_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);

% بررسی خروجی
if exitflag == 1 || exitflag == 0 % 0 یعنی به حداکثر تکرار رسیده اما جواب قابل قبول است
    u = z_opt(1:n_u);
    delta_clf = z_opt(n_z);
else
    % استراتژی پشتیبان (Fallback) اگر QP شکست خورد
    % استفاده از کنترلر PD ساده برای حفظ پایداری عددی
    u_pd = -(20*e + 5*de);
    u = max(min(u_pd, u_max), -u_max);
    delta_clf = 0;
    % نمایش هشدار در کنسول برای عیب‌یابی (اختیاری)
    % fprintf('Warning: QP failed at t=%.3f, Exitflag=%d\n', t, exitflag);
end

end

% ---------------------------------------------------------------------
% توابع کمکی
% ---------------------------------------------------------------------

function [y, J_y, Jdot_y_dq] = get_ActualOutputs(q, dq)
    % خروجی‌ها: زوایای ۴ مفصل ران و زانو
    % ترتیب q در ربات Rabbit: [x, y, theta, q1, q2, q3, q4]
    H_0 = [zeros(4, 3), eye(4)];
    y = H_0 * q;
    J_y = H_0;
    Jdot_y_dq = zeros(4, 1); 
end

function [y_d, dy_d, ddy_d] = get_DesiredOutputs(t, x0, p)
    % GET_DESIREDOUTPUTS  Generates trajectory consistent with x0
    % t  : Current time
    % x0 : The 14x1 initial state vector from make_initial_state
    % p  : Parameter vector
    
    % استخراج پیکربندی اولیه از x0
    q0 = x0(1:7);
    dq0 = x0(8:14);
    
    % مفصل‌های مورد نظر برای ردیابی (q1, q2, q3, q4)
    % طبق make_initial_state، این‌ها اندیس‌های 4 تا 7 هستند
    y_init_vals = q0(4:7); 
    
    % مشتقات اولیه برای شروع نرم (بسیار مهم!)
    % از آنجایی که در make_initial_state سرعت‌ها صفر نیستند (به دلیل اصلاح J*dq=0)
    % ما باید سرعت اولیه را از x0 بگیریم
    dy_init_vals = dq0(4:7); 
    
    % پارامترهای نوسان (دامنه را بسته به نیاز تنظیم کنید)
    Amp = [0.1; 0.15; 0.1; 0.15]; 
    omega = 2 * pi * 1.5; 
    
    % --- فرمول اصلاح شده با در نظر گرفتن سرعت اولیه ---
    % برای اینکه در t=0 هم موقعیت و هم سرعت با x0 یکی باشد:
    % y_d(t) = y_init + (dy_init/omega) * sin(omega*t) + Amp * (1 - cos(omega*t))
    
    % سینوس برای حفظ سرعت اولیه (dy_init)
    term1 = (dy_init_vals ./ omega) .* sin(omega * t);
    
    % (1-cos) برای شروع نرم از موقعیت y_init
    term2 = Amp .* (1 - cos(omega * t));
    
    y_d   = y_init_vals + term1 + term2;
    dy_d  = dy_init_vals .* cos(omega * t) + (Amp .* omega) .* sin(omega * t);
    ddy_d = -(dy_init_vals .* omega) .* sin(omega * t) + (Amp .* omega.^2) .* cos(omega * t);
end


