function coeffs = make_coeffs(x, p)
%MAKE_COEFFS Initial B-spline control points for the virtual constraints.
%
% Smooth interpolation from x's actuated joints to a periodicity-motivated
% touchdown target: for a periodic step the pose right before impact
% (pre-relabel) should look like the current pose with stance/swing roles
% swapped -- exactly what relabel_state computes.
%
% Because a clamped B-spline's first control point IS its value at s=0,
% this also guarantees spline(s=0) == x's actuated joints exactly. Rebuild
% the coefficients whenever x changes so that invariant is preserved.

    x_relabeled = relabel_state(x);
    q_end_guess = x_relabeled(4:7);

    coeffs = zeros(p.nu, p.n_coeffs);
    for i = 1:p.nu
        q_start_i = x(i+3);          % act_idx = 4:7, offset by 3
        q_end_i   = q_end_guess(i);
        coeffs(i,:) = linspace(q_start_i, q_end_i, p.n_coeffs);
    end
end
