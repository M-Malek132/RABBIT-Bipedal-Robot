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

% WHICH ROBOT. [] is today's 74 kg model (M.m / V.m / G.m). A scalar lambda in
% [0, 1] blends in the 30 kg model this repository used before 2026-09-02:
% 0 is that robot, 1 is today's, and anything between has linearly blended
% masses and inertias -- exact, since the dynamics are linear in them (ch3_mvg).
% Only ch3_model_homotopy sets it.
p.model_blend = [];

%% phase variable theta
p.c_theta = [0 0 1 1 0.5 0 0];   % theta = qt + q1 + q2/2 = c*q
p.theta_minus = -0.15;
p.theta_plus  =  0.30;

% THE PHASE ENDPOINTS AS DECISION VARIABLES. With both legs of 0.5 m links,
% theta = qt + q1 + q2/2 is EXACTLY the direction of the hip-to-foot line
% (Dynamics/DOCS.md), and at the strike both feet are on the ground with the
% stance leg at theta_plus and the landing leg at theta_minus (periodicity).
% So the stride is fixed by the hip height at impact alone:
%
%       L_step = h_hip(t_N) * (tan theta_plus - tan theta_minus)
%
% which every stored gait satisfies to ~1e-7. With theta_minus, theta_plus
% held constant the only way to shorten the stride is to lower the hip at
% impact, and the hip band (row 8) floors that: at 0.925 - 0.045 m the stride
% cannot go below 0.405 m. free_theta = true appends [theta_minus; theta_plus]
% to the collocation vector (ch3_col_pack) inside theta_bounds, so the stride
% becomes a design variable. OFF by default: nothing solved so far changes.
p.free_theta   = false;
p.theta_bounds = [-0.40 -0.02; ...  % theta_minus range [rad]
                   0.10  0.50];      % theta_plus  range [rad]

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
p.limits.clearance_max = 0.15;  % swing-foot height CEILING, whole step [m]
% CEILING on step length, gated OFF (p.limits.enable.step_len_max) so nothing
% solved so far changes behaviour. Turn it on to ask for a SHORT stride, which
% the objective will not give on its own: torque-squared per distance makes a
% longer stride cheaper, so every gait solved on this model sits at
% L = 0.438 m. A short stride is what a slow gait needs -- on a long one the
% contact runs out, Fz sags to its floor and the friction demand pins at mu_s,
% which is what walls the speed march at 1.15 m/s (row 19,
% ch3_col_constraints).
p.limits.step_len_max  = 0.30;  % L_step <= this, when enabled       [m]
% Objective: torque-squared per unit DISTANCE (true) or per STEP (false). The
% division by step length is what keeps the stride from collapsing, but it also
% opposes the ceiling above -- with both stride rows on, the constraints box the
% stride and the normalization only gets in the way. Leave TRUE unless marching
% the stride down, and never turn it off with the ceiling off.
p.cost_normalize = true;

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
                         'clearance_max', true, ...
                         'height',  true, ...
                         'swing_clear', true, ...  % NIC3
                         'liftoff',     true, ...  % NEC2
                         'impact',      true, ...  % NEC3
                         'hzd',         true, ...  % NEC4 + NEC5
                         'phase_mono',  true, ...  % HH6
                         'decoupling',  true, ...  % HH2
                         'step_len_max', false);   % OFF: see step_len_max below

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

p.scale_problem = false;        % fmincon ScaleProblem. FALSE keeps the reported
                                % Feasibility directly comparable to
                                % ConstraintTolerance (see ch3_col_solve); set
                                % TRUE when a march stalls at exitflag -2 on a
                                % badly scaled rung.

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

% END A STEP AT ITS FIRST CONTACT-INVALID SAMPLE (sampled control only). The
% simulator pins the stance foot, so a step whose true contact force lifts off
% (Fz <= 0) or slips (|Fx| > mu_s Fz) still "completes"; ch3_validity scores
% that afterwards. With this on, the step stops at that sample instead and is
% reported as failed, so the steps a run completes ARE its valid steps and no
% metric is taken over a trajectory the ground could not have produced. OFF by
% default so every stored table reproduces.
p.stop_on_invalid = false;

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
