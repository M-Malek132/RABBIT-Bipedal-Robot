function p = ch3_params(varargin)
%CH3_PARAMS  Single source of truth for the Chapter-3 pipeline.
%
%   p = ch3_params()                 default parameter struct
%   p = ch3_params('name',value,...) defaults with overrides

p = struct();

%% model
p.nq   = 7;                     % generalized coordinates
p.nu   = 4;                     % actuators
p.nx   = 2*p.nq;                % full state dimension (14)
p.iact = (4:7)';                % actuated coordinate indices
p.ny   = numel(p.iact);         % number of outputs (4)
p.g0   = 9.8062;                % gravity used in the symbolic derivation
p.mass = 74;                    % total mass [kg]

p.H = zeros(p.ny, p.nq);
p.H(:, p.iact) = eye(p.ny);

%% phase variable theta
p.c_theta = [0 0 1 1 0.5 0 0];   % theta = qt + q1 + q2/2 = c*q
p.theta_minus = -0.15;
p.theta_plus  =  0.30;

%% virtual constraint parametrization
p.basis   = 'bspline';          % 'bezier' | 'bspline'
p.bez_deg = 5;                  % Bezier degree M
p.n_ctrl  = p.bez_deg + 1;      % alpha columns per output (M+1)
p.bsp_deg = 3;                  % degree used only when basis = 'bspline'

%% controller
p.controller = 'clfqp_con';     % 'iolin_pd' | 'clfqp' | 'clfqp_con' | 'ff'

p.eps = 0.50;
p.Kp  = eye(p.ny) * 100;
p.Kd  = eye(p.ny) *  20;

%% RES-CLF
p.clf_construction = 'lyap';    % 'care' | 'lyap'
p.Q_clf = eye(2*p.ny);          % Q in the (C)ARE / Lyapunov equation
p.clf_slack_penalty = 1e6;      % "p" multiplying delta^2 in stage 8

%% gait / speed targets
p.v_des        = 1.2;           % NEC1 average walking rate [m/s]
p.step_len_min = 0.15;          % floor on step length, kills "step in place"
p.qt_range     = [-1.0 1.0];    % [qt_min qt_max] torso pitch box [rad]
p.enforce_nec1 = true;          % gate on the NEC1 speed equality

%% Table 3.1 physical realizability limits
p.limits = struct();
p.limits.u_max       = 120;     % |u_i| <= u_max                    [Nm]
p.limits.impulse_max = 15;      % ||impact impulse||_2 <= impulse_max [Ns]
p.limits.mu_s        = 0.4;     % |Fx| <= mu_s * Fz     NIC2        [-]
p.limits.Fz_min      = 50;      % Fz >= Fz_min          NIC1        [N]
p.limits.clearance   = 0.05;    % swing-foot height at mid-step     [m]

%% Section 6.3.4 constraint set (NIC / NEC)
p.limits.sw_clear_min   = 1e-3; % strict swing-foot clearance, interior  [m]
p.limits.sw_strike_rate = 0.05; % foot must be descending at strike    [m/s]
p.limits.liftoff_rate   = 0.01; % d/dt(new swing-foot height) at t=0+  [m/s]
p.limits.impulse_z_min  = 0;    % Iz >= this (compressive)              [Ns]
p.limits.mu_s_impact    = [];   % |Ix| <= this * Iz; [] -> p.limits.mu_s   [-]
p.limits.thetadot_min   = 0.10; % thetadot >= this over the step      [rad/s]
p.limits.dec_min        = 1e-3; % sigma_min(LgLf y) >= this
p.limits.hip_h       = 0.90;    % centre of the hip-height band      [m]
p.limits.hip_h_tol   = 0.02;    % half-width (so 4 cm of bob total)  [m]

p.limits.enable = struct('torque',  true, ...
                         'impulse', true, ...
                         'friction',true, ...   % NIC2
                         'grf',     true, ...   % NIC1
                         'clearance', true, ...
                         'height',  true, ...
                         'swing_clear', true, ...  % NIC3
                         'liftoff',     true, ...  % NEC2
                         'impact',      true, ...  % NEC3
                         'hzd',         true, ...  % NEC4 + NEC5
                         'phase_mono',  true, ...  % HH6
                         'decoupling',  true);     % HH2

%% hybrid zero dynamics (NEC4/5)
p.hzd_grid       = 161;         % reporting quadrature grid
p.hzd_grid_solve = 41;          % constraint quadrature grid

p.hzd_tol = struct('delta_margin', 1e-3, ...   % delta^2 <= 1 - this
                   'delta_min',    1e-8, ...   % delta^2 >= this
                   'zeta_margin',  1e-6);      % zeta* clears V_max/delta^2 by this

%% direct collocation
p.N_nodes  = 41;                % Hermite-Simpson nodes per step
p.T_min    = 0.20;              % step duration bounds [s]
p.T_max    = 1.50;
p.dq_max   = 20;                % box on joint velocities

p.max_iter      = 300;
p.max_fun_evals = 3e5;

p.verify_tol = 1e-3;            % ch3_col_verify rejection tolerance

p.checkpoint_file  = '';        % path to resume a solve; '' disables
p.checkpoint_every = 10;

p.seed_knee_lift = 1.2;         % mid-step swing-knee flexion in seed [rad]

%% integration
p.ode_reltol = 1e-8;
p.ode_abstol = 1e-9;
p.ode_maxstep = 5e-3;
p.guard_min_time = 0.05;        % ignore guard crossings before this [s]

p.control_dt = 0;               % 0 = evaluate feedback continuously

%% overrides
for k = 1:2:numel(varargin)
    field = varargin{k};
    if ~isfield(p, field)
        error('ch3_params:unknownField', 'Unknown parameter "%s".', field);
    end
    p.(field) = varargin{k+1};
end

% Keep derived quantities consistent if the caller changed bez_deg.
p.n_ctrl = p.bez_deg + 1;
p.ny     = numel(p.iact);

end
