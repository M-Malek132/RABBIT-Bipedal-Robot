function u = rabbit_pd_controller(x,p)

q  = x(1:7);
dq = x(8:14);

% actuated joints
qa  = q(4:7);
dqa = dq(4:7);

% desired configuration
qa_des  = p.qa_des;
dqa_des = zeros(4,1);

% gains
Kp = diag([100 100 80 80]);
Kd = diag([20 20 15 15]);

% PD control
u = -Kp*(qa - qa_des) ...
    -Kd*(dqa - dqa_des);

end
