function [y,dy] = rabbit_virtual_constraints(q,dq,p)

theta  = q(3);
dtheta = dq(3);

s = (theta - p.theta0) / (p.thetaf - p.theta0);

hd  = bezier(s,p.alpha);
dhd = bezier_derivative(s,p.alpha) ...
      * dtheta/(p.thetaf - p.theta0);

qa  = q(4:7);
dqa = dq(4:7);

y  = qa  - hd;
dy = dqa - dhd;

end