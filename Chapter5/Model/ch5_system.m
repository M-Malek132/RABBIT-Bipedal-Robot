function p = ch5_system(p, explicit)
%CH5_SYSTEM  Merge the selected plant's block into the parameter struct.
%
%   p = ch5_system(p)
%   p = ch5_system(p, explicit)   explicit = cellstr of names the caller set
%
% Chapter 5 validates one method on two deliberately different plants, and the
% differences are the point:
%
%                      springmass (Fig 5.1)        pendulum (Fig 5.2)
%   ------------------------------------------------------------------------
%   dynamics           LINEAR                      NONLINEAR
%   states  nx         6                           8
%   inputs  nu         1                           2
%   outputs ny         1                           2
%   relative degree    6                           4
%   DOF / underact.    3 / 2                       4 / 2
%   constraint         x3 <= x3max  (position)     py2 >= p2min (end effector)
%
% Both have relative degree EQUAL TO nx/ny, so both are fully input-output
% linearizable with NO ZERO DYNAMICS. That is a convenience of the validation
% systems, not of the method: eta is the whole state, so nothing can hide in an
% internal mode and the plots show everything there is.
%
% ONE STRUCTURAL FACT WORTH STATING UP FRONT, because it is why these two
% plants were chosen and why the barrier is hard. In both systems the
% CONSTRAINED QUANTITY IS THE CONTROLLED QUANTITY -- x3 is both the tracking
% output and the thing bounded; the end-effector height is a function of the
% same theta the outputs track. So the safety constraint inherits the FULL
% relative degree of the plant. There is no cheap reformulation in which the
% barrier is relative degree 1, which is exactly the situation Section 5.1's
% reciprocal CBF cannot address and Section 5.2 exists to fix.
%
% Defaults set here NEVER overwrite a name the caller passed to ch5_params;
% that is what `explicit` is for.
%
% See also CH5_PARAMS, CH5_SPRINGMASS, CH5_PENDULUM.

if nargin < 2, explicit = {}; end

switch lower(p.system)

    %% =================================================== serial spring mass
    case 'springmass'

        % Eq (5.32)-(5.34). The thesis does not publish m_i or k, so these are
        % chosen -- unit masses and unit stiffness -- and they are the reason
        % the reproduced figure matches Fig. 5.3 in magnitude:
        %
        %   * the chain's natural frequencies are 0.445, 1.25, 1.80 rad/s, so a
        %     30 s window holds the whole transient, as Fig. 5.3 shows;
        %   * L_g L_f^5 y = k^2/(m1 m2 m3) = 1, so the decoupling "matrix" is
        %     perfectly conditioned and the plotted forces are the controller's
        %     doing rather than an inversion artifact;
        %   * u -> 0 at the target, because x1 = x2 = x3 relaxes every spring.
        d = struct();
        d.m = [1 1 1];              % cart masses            [kg]
        d.k = 1;                    % spring stiffness       [N/m]

        d.x0  = zeros(6,1);         % all carts at the origin, at rest
        d.x3d = 3;                  % desired position of cart 3     [m]

        % Fig. 5.3 sweeps the constraint while HOLDING THE POLES FIXED, which
        % is the figure's actual claim: the pole locations encode the response,
        % the constraint value encodes the limit, and the two are separable.
        %   (b) x3max - x3d = 15 cm
        %   (c) x3max - x3d =  0 cm   <- target sits exactly on the boundary
        d.constraint = struct('type', 'x3_max', 'value', 3.15);

        d.T = 30;

        % pb = 0.12 * [10 11 12 13 14 15], straight from Section 5.2.3.
        d.ecbf_poles = 0.12 * (10:15);

        % CLF RATE. The thesis does not publish lambda, so it is set here
        % against a number the thesis DOES publish: Section 5.2.3 reports that
        % the CLF-QP baseline overshoots to max(x3) = 3.27 m. Measured on this
        % model, lambda = 0.20 gives 3.29 -- the closest of the values swept
        % (0.18 -> 3.36, 0.22 -> 3.22), and it puts the whole transient inside
        % the 30 s window with peak forces of order 1 N as Fig. 5.3 shows.
        %
        % The rate is not free above this. The min-norm CLF-QP applies
        % ||mu|| = psi/||LgV||, so demanding a rate the geometry cannot supply
        % divides by a small number: lambda = 0.3 already costs 24 N and
        % lambda = 1.0 diverges. That ceiling is a property of the min-norm
        % formulation, not of the plant, and it is why lambda is a per-system
        % default rather than one shared constant.
        d.clf_lambda = 0.20;
        d.clf_Q      = eye(6);

        d.sys = struct( ...
            'name',   'springmass', ...
            'pretty', 'Serial spring-mass (relative degree 6)', ...
            'nx', 6, 'nu', 1, 'ny', 1, 'r', 6, 'rb', 6, ...
            'nq', 3, ...
            'x_labels', {{'x_1','x_2','x_3','$\dot x_1$','$\dot x_2$','$\dot x_3$'}}, ...
            'u_labels', {{'u'}}, ...
            'u_unit',   'N', ...
            'y_labels', {{'x_3 - x_{3d}'}});

    %% ============================ two-link pendulum with elastic actuators
    case 'pendulum'

        % Eq (5.36)-(5.37). The thesis publishes only "the links have unit
        % length"; every other number below is a choice, recorded here rather
        % than buried, so the reproduced Fig. 5.4 is honest about what it is.
        d = struct();
        d.m  = [1 1];               % link masses                   [kg]
        d.l  = [1 1];               % link lengths                  [m]   (given)
        d.lc = [0.5 0.5];           % COM along each link           [m]
        d.I  = [1 1]/12;            % link inertia about COM (uniform rod)
        d.g  = 9.81;

        % MOTOR / SERIES-ELASTIC BLOCK. These are what make the system
        % relative degree 4 instead of 2: tau reaches the link only through the
        % spring, i.e. through TWO extra integrators per joint.
        %
        % k is a compliance, and the softer it is the further the actuator sits
        % from the load. 100 Nm/rad with Jm = 0.1 puts the motor-spring mode at
        % sqrt(k/Jm) = 31.6 rad/s, comfortably above the arm's own ~3 rad/s and
        % comfortably below the 1 kHz control rate -- a real SEA regime rather
        % than a stiff-limit one that would collapse back to relative degree 2.
        d.Jm = 0.1;                 % motor inertia                 [kg m^2]
        d.k  = 100;                 % motor/joint spring stiffness  [Nm/rad]
        d.xi = 1.0;                 % joint damping of (5.37)       [Nm s/rad]

        % Fig. 5.4: theta1 from -pi to +pi, theta2 held at 0.
        %
        % READ THAT AGAIN -- both are the arm pointing STRAIGHT UP, and they
        % differ by a full turn. The task is to rotate link 1 all the way
        % around, and the only way through is theta1 = 0, where an unfolded arm
        % puts the end effector at py = -2 m. The safety constraint py >= p2min
        % is therefore not a mild bound near the trajectory; it is
        % incompatible with theta2 = 0 over part of the maneuver and FORCES the
        % arm to fold against what the CLF is asking for. That conflict is what
        % Fig. 5.4's "aggressively move the links" is describing.
        d.theta0 = [-pi; 0];
        d.thetad = [ pi; 0];

        %   (b) p2min = -1.0 m
        %   (c) p2min = -0.5 m   <- tighter, and visibly more expensive
        d.constraint = struct('type', 'py_min', 'value', -1.0);

        % 15 s, the window Fig. 5.4 spans.
        %
        % The ECBF runs do not always finish the maneuver inside it, and that
        % is the result rather than a reason to lengthen the window. The
        % barrier row is hard and the CLF row is slacked, so where the two
        % conflict -- and here they conflict for seconds, while theta2 folds
        % and unfolds -- tracking is DEFERRED, not abandoned. How far it gets
        % by t = 15 is not a robust number anyway: once the barrier is active
        % the closed loop is a discrete-time sliding mode and the trajectory is
        % integrator-dependent (see the note on d.integrator). What is robust,
        % and what the chapter claims, is that the safe set is never left.
        d.T = 15;

        % pb = [5 5.5 6 10], from Section 5.2.3.
        d.ecbf_poles = [5 5.5 6 10];

        % CLF WEIGHT AND RATE, both chosen by measurement rather than taste.
        %
        % Q = 1000 I is not cosmetic. The task is a 2*pi rotation, so y1(0) =
        % -2*pi, and a min-norm CLF-QP riding Vdot = -lambda V follows whatever
        % path the level sets of V dictate. With Q = I the CARE puts the
        % closed-loop poles on the unit circle and the arm does not finish
        % within 15 s -- swept: Q = 10 leaves |theta1 - pi| = 0.86, Q = 100
        % leaves 0.12, Q = 1000 leaves 0.020. Q = 3000 buys nothing further.
        %
        % lambda = 1.0 with that Q completes the maneuver with peak torques
        % near 100 Nm, matching the scale drawn in Fig. 5.4a. Faster is worse
        % in the same way as on the springmass: lambda = 2 costs 280 Nm.
        d.clf_lambda = 1.0;
        d.clf_Q      = 1000 * eye(8);

        % 2 kHz. 1 kHz already tracks the barrier well, but sampled-data
        % enforcement leaves an O(dt) excursion past the boundary -- the ECBF
        % row is imposed at sample instants and the plant moves between them.
        % Measured at p2min = -1.0: dt = 2 ms gives h_min = -5.1e-3, 1 ms gives
        % -3.7e-4, and 0.5 ms gives +6.4e-4, i.e. STRICTLY SAFE. Since the
        % chapter's claim is forward invariance, the default is the first rate
        % at which the reproduction actually exhibits it rather than nearly
        % does. ch5_test_ecbf measures this excursion rather than assuming it
        % away.
        d.control_dt = 5e-4;

        % Fixed-step RK4 rather than the ode45 default, and the reason is not
        % speed alone.
        %
        % Once the barrier row is active this closed loop is a discrete-time
        % sliding mode, and an INDIVIDUAL TRAJECTORY is not reproducible across
        % integrators: rk4 and ode45 separate by ~0.5 rad within 3 s (measured
        % in ch5_test_ecbf). Neither is "the" answer, and paying ode45's
        % adaptive-step cost buys no fidelity in a quantity that is not
        % reproducible anyway. What IS reproducible -- and is what the chapter
        % claims -- is that the safe set is not left; both integrators agree on
        % that, and ch5_test_ecbf asserts it for both.
        %
        % At 10 substeps of the 0.5 ms period, RK4 takes 50 us steps against a
        % fastest plant mode of 31.6 rad/s: about 4000 steps per period.
        d.integrator = 'rk4';

        d.sys = struct( ...
            'name',   'pendulum', ...
            'pretty', '2-link pendulum, elastic actuators (relative degree 4)', ...
            'nx', 8, 'nu', 2, 'ny', 2, 'r', 4, 'rb', 4, ...
            'nq', 4, ...
            'x_labels', {{'\theta_1','\theta_2','\theta^m_1','\theta^m_2', ...
                          '$\dot\theta_1$','$\dot\theta_2$', ...
                          '$\dot\theta^m_1$','$\dot\theta^m_2$'}}, ...
            'u_labels', {{'\tau_1','\tau_2'}}, ...
            'u_unit',   'Nm', ...
            'y_labels', {{'y_1 = \theta_1 - \theta_{1d}', ...
                          'y_2 = \theta_2 - \theta_{2d}'}});

    otherwise
        error('ch5_system:unknownSystem', ...
              'Unknown p.system "%s" (expected springmass|pendulum).', p.system);
end

%% ---- merge, without stepping on anything the caller asked for explicitly
% p.sys and the plant block go under their own names so the two plants never
% collide; the shared names (T, control_dt, ecbf poles) are filled in only if
% the caller left them alone.
shared  = {'sys','T','control_dt','ecbf_poles','clf_lambda','clf_Q','integrator'};
p.plant = rmfield(d, intersect(fieldnames(d), shared));
p.sys   = d.sys;

if ~was_set(explicit, 'T') || isempty(p.T)
    p.T = d.T;
end
if isfield(d, 'control_dt') && ~was_set(explicit, 'control_dt')
    p.control_dt = d.control_dt;
end
if isfield(d, 'integrator') && ~was_set(explicit, 'integrator')
    p.integrator = d.integrator;
end
if isempty(p.ecbf.poles) && ~was_set(explicit, 'ecbf.poles')
    p.ecbf.poles = d.ecbf_poles;
end
if isempty(p.clf.lambda) && ~was_set(explicit, 'clf.lambda')
    p.clf.lambda = d.clf_lambda;
end
if isempty(p.clf.Q)
    p.clf.Q = d.clf_Q;
end

% Pack the pendulum parameters once, here, rather than in every ODE right-hand
% side evaluation. The generated Lie-derivative code takes this vector.
if strcmpi(p.system, 'pendulum')
    p.plant.pv = ch5_pend_pv(p);
end

% p.constraint overrides the plant default: a scalar swaps the level only, a
% struct swaps the whole spec. See the note in ch5_params.
if ~isempty(p.constraint)
    if isnumeric(p.constraint) && isscalar(p.constraint)
        p.plant.constraint.value = p.constraint;
    elseif isstruct(p.constraint)
        p.plant.constraint = p.constraint;
    else
        error('ch5_system:constraint', ...
              'p.constraint must be [] , a scalar level, or a spec struct.');
    end
end
p.constraint = p.plant.constraint;      % resolved, so callers can read one place

end

% ---------------------------------------------------------------------------
function tf = was_set(explicit, name)
tf = any(strcmpi(explicit, name));
end
