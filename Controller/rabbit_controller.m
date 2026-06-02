function u = rabbit_controller(t, x, params)

nq = 7;

q  = x(1:nq);
dq = x(nq+1:end);

theta = q(4) + 0.0*q(6);

theta0 = -0.3;
thetaf = 0.3;

s = (theta - theta0)/(thetaf - theta0);
s = min(max(s,0),1);

hd = desired_gait(s);

y  = q(4:7) - hd;
dy = dq(4:7);

Kp = diag([200 200 150 150]);
Kd = diag([30 30 20 20]);

u = -Kp*y - Kd*dy;

end
