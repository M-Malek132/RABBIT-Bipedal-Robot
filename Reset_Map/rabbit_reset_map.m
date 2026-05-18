function x_plus = rabbit_reset_map(x_minus, param)

q  = x_minus(1:7);
dq = x_minus(8:14);

% --- kinematics ---
[~, p_sw, ~, ~, ~, ~] = rabbit_kinematics(q,packParameters(param));

% --- vertical translation ---
q(1) = q(1) - p_sw(1);
q(2) = q(2) - p_sw(2);


% --- compute impact velocity ---
x_tmp = rabbit_impact_map([q; dq], packParameters(param));
dq_plus = x_tmp(8:14);

% --- leg relabeling ---
q_plus = q;
q_plus([4 5 6 7]) = q([6 7 4 5]);

dq_tmp = dq_plus;
dq_plus([4 5 6 7]) = dq_tmp([6 7 4 5]);

% --- ensure column vectors ---
q_plus  = q_plus(:);
dq_plus = dq_plus(:);

% --- assemble state ---
x_plus = [q_plus; dq_plus];

end
