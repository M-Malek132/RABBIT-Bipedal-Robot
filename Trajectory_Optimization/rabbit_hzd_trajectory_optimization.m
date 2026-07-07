%% =========================================================
% HZD PARAMETERS
%% =========================================================
n_coeffs = 6; % Number of control points per joint (degree 3-5)
nu = 4;       % Actuated joints
z0 = randn(nu * n_coeffs, 1); % Initial guess for coefficients

%% =========================================================
% OPTIMIZATION
%% =========================================================
options = optimoptions('fmincon', 'Display','iter', 'Algorithm','sqp');

[z_opt, fval] = fmincon(...
    @(z) hzd_cost(z, p), ...      % Cost: Minimize energy/effort
    z0, [], [], [], [], [], [], ...
    @(z) hzd_constraints(z, p), ... % Constraints: Periodicity & Stability
    options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HZD COST FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = hzd_cost(z, p)
    % 1. Unpack coefficients into joint trajectories
    % 2. Simulate one full gait cycle (using ode45)
    % 3. J = Integral of (Torque^2) over one cycle
    
    coeffs = reshape(z, [p.nu, p.n_coeffs]);
    [~, ~, total_torque_sq] = simulate_hzd_gait(coeffs, p);
    J = total_torque_sq;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HZD CONSTRAINTS (The "Meat" of HZD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c, ceq] = hzd_constraints(z, p)
    coeffs = reshape(z, [p.nu, p.n_coeffs]);
    
    % 1. SIMULATION
    % We simulate the robot using a PD controller tracking the B-spline
    [x_start, x_end] = simulate_hzd_gait(coeffs, p);
    
    % 2. PERIODICITY (The Poincaré Map)
    % After impact (relabeling), the state must match the start
    x_next = rabbit_impact_map(x_end);
    
    % This is the primary equality constraint
    ceq = x_next - x_start; 
    
    % 3. STABILITY (Optional but recommended)
    % To ensure the gait isn't just periodic but also stable,
    % we want the eigenvalues of the Monodromy matrix to be < 1.
    % Since fmincon handles equality/inequality well, we can add:
    % c = max(0, eigenvalues(M) - 0.95); 
    
    c = []; 
end
