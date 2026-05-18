function u = rabbit_fl_controller(x,p)

q  = x(1:7);
dq = x(8:14);

[Dc,Hc,Bc] = rabbit_constrained_dynamics(q,dq,p);

[y,dy,Hq,dHq] = rabbit_outputs(q,dq,p);

A = Hq*(Dc\Bc);

b = -Hq*(Dc\Hc) + dHq*dq;

Kp = 50*eye(4);
Kd = 15*eye(4);

v = -Kp*y - Kd*dy;

u = A\(v - b);

end
