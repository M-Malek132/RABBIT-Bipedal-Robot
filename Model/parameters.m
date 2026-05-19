function p = parameters()
% PARAMETERS  Return structure of physical parameters for RABBIT robot.

% Masses [kg]
p.mT = 10;     % torso
p.m1 = 5;      % thigh
p.m2 = 5;      % shank

% Link lengths [m]
p.l1 = 0.5;
p.l2 = 0.5;
p.lt = 0.75;

% Inertias [kg·m²]
p.I1 = 0.1;
p.I2 = 0.1;
p.IT = 0.2;

% Gravity [m/s²]
p.g  = 9.81;

% --- Generate Random Discrete Stepping Stones ---
num_stones = 14;
stones = zeros(num_stones + 2, 2);

% Initial ground (covers starting stance foot safely)
stones(1, :) = [-0.5, 0.5];
current_x = 0.5;

% Randomly generate intermediate stones
for i = 1:num_stones
    % Define random gap and width bounds
    gap = 0.05 + 0.15 * rand();   % Random gap between 0.05m and 0.20m
    width = 0.3 + 0.3 * rand();   % Random width between 0.30m and 0.60m
    
    start_x = current_x + gap;
    end_x = start_x + width;
    
    stones(i + 1, :) = [start_x, end_x];
    current_x = end_x;
end

% Solid ground for the rest of the simulation
final_gap = 0.05 + 0.15 * rand();
stones(end, :) = [current_x + final_gap, current_x + final_gap + 6.0];

p.stones = stones;

end
