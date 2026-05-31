function [hd, dhd_ds] = desired_gait_bspline(s, ControlPoints, n, p)
% Evaluate desired joint trajectories using B-spline
%
% Inputs:
%   s:             normalized phase variable [0,1]
%   ControlPoints: (n+1) x 4 matrix of B-spline control points
%                  Column j = control points for joint j
%   n:             number of data points (control points = n+1)
%   p:             B-spline degree
%
% Outputs:
%   hd:     4x1 desired joint angles [q1;q2;q3;q4]
%   dhd_ds: 4x1 derivative w.r.t phase (computed numerically)

    % Clamp s to valid range
    s = min(max(s, 0), 1);
    
    % Evaluate B-spline basis at s
    N = BSpline(n, p, s);  % 1 x (n+1) row vector
    
    % Compute desired angles for all 4 joints
    hd = zeros(4, 1);
    for j = 1:4
        hd(j) = N * ControlPoints(:, j);
    end
    
    % Compute derivative numerically (central difference)
    if nargout > 1
        ds = 1e-4;  % small perturbation for numerical derivative
        
        % Evaluate at s+ds and s-ds (handle boundaries)
        s_plus  = min(s + ds, 1);
        s_minus = max(s - ds, 0);
        
        N_plus  = BSpline(n, p, s_plus);
        N_minus = BSpline(n, p, s_minus);
        
        dhd_ds = zeros(4, 1);
        for j = 1:4
            h_plus  = N_plus  * ControlPoints(:, j);
            h_minus = N_minus * ControlPoints(:, j);
            dhd_ds(j) = (h_plus - h_minus) / (s_plus - s_minus);
        end
    end
end