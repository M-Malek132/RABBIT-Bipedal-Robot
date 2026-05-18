function hd = desired_gait(s)

hd = zeros(4,1);

% stance leg
hd(1) = -0.3 + 0.4*s;
hd(2) = 0.6 - 0.25*sin(pi*s);

% swing leg
hd(3) = -0.8 + 1.2*s;
hd(4) = 0.5 + 0.7*sin(pi*s);

end
