function u = rabbit_hzd_controller(x,p)

q  = x(1:7);
dq = x(8:14);

theta  = q(3);
dtheta = dq(3);

s = (theta - p.theta0)/(p.thetaf - p.theta0);

hd  = bezier(s,p.alpha);
dhd = bezier_derivative(s,p.alpha)*dtheta/(p.thetaf - p.theta0);

qa  = q(4:7);
dqa = dq(4:7);

y  = qa - hd;
dy = dqa - dhd;

[Dc,Hc,Bc] = rabbit_constrained_dynamics(q,dq,p);

Hq = [zeros(4,3) eye(4)];

A = Hq*(Dc\Bc);
b = -Hq*(Dc\Hc);

Kp = 80*eye(4);
Kd = 20*eye(4);

v = -Kp*y - Kd*dy;

u = A\(v - b);

end
