function dx = col_dynamics(x, u, lam, p)
%COL_DYNAMICS  Single-support state derivative with EXPLICIT contact force.
%
%   dx = col_dynamics(x, u, lam, p)     x = [q; dq], dx = [dq; ddq]
%
% Unlike rabbit_constrained_dynamics (which SOLVES the KKT system for the
% contact force), the collocation transcription carries the stance-foot force
% lam = [Fx; Fz] as an explicit decision variable and evaluates the plain
% equation of motion
%
%       M(q) ddq + C + G = B u + J_st(q)' lam.
%
% The holonomic contact constraint (foot pinned) is then imposed SEPARATELY as
% path constraints in col_constraints (position at every node, velocity at the
% first node), which is what keeps the stance foot from drifting -- the flaw
% the shooting simulator documents. This split is the standard index-reduced
% form for constrained-multibody direct collocation.

    nq = p.nq;
    q  = x(1:nq);
    v  = x(nq+1:2*nq);

    Mm = M(q);
    Cv = V([q; v]);
    Gv = G(q);
    B  = input_matrix();
    Js = J_st(q);

    ddq = Mm \ (B*u + Js.'*lam - Cv - Gv);
    dx  = [v; ddq];
end
