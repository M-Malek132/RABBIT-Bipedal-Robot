function [x0] = make_initial_state()
% MAKE_INITIAL_STATE Construct a consistent initial state for RABBIT.
% No longer requires parameter structs.

    nq = 7;
    x0 = zeros(2 * nq, 1);

    %% 1. Define initial configuration
    % Using arbitrary constants as you have defined your own kinematics
    qt = 0.1; q1 = -0.3; q2 = 0.6; q3 = -1.0; q4 = 0.6;
    
    % We need an initial (px, pz). 
    % If your kinematics are relative, you might set px=0, pz=0 to start.
    px = 0; pz = 0; 
    
    q0 = [px; pz; qt; q1; q2; q3; q4];
    x0(1:nq) = q0;

    %% 2. Initial velocities
    dq0 = zeros(nq, 1);
    dq0(3) = 0.3; % Torso velocity

    %% 3. Enforce stance foot contact: J_stance * dq = 0
    % Assuming J_stance(q) exists independently of params
    J = J_st(q0); 
    
    % Project velocity onto the null space of the constraint Jacobian
    dq0_corrected = (eye(nq) - pinv(J) * J) * dq0;
    x0(nq + 1:end) = dq0_corrected;
    
    fprintf('Initial state initialized. Foot constraint enforced.\n');
end
