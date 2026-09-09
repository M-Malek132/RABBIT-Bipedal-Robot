function [Mq, Vv, Gv, tag] = ch4_case_dynamics(s, q, dq)
%CH4_CASE_DYNAMICS  Dispatch to the independently re-derived per-case M/V/G.
%
%   [Mq, Vv, Gv, tag] = ch4_case_dynamics(s, q, dq)
%   Mq                 = ch4_case_dynamics(s, q)          (V not needed)
%
% Chapter 4's mass_scale perturbation is no longer applied as s*M(q) inline
% inside ch4_control_affine/ch4_impact. Each scale the chapter actually uses
% has its own M_<tag>/V_<tag>/G_<tag>, rerun from scratch through the
% Lagrangian trace derivation by rabbit_generate_case_dynamics (see that
% file's header for why: it does not trust the "uncertainty is linear in
% mass, so M -> sM" argument on faith, it reruns the derivation and lets
% ch4_test_model compare the two numerically).
%
% s = 1 is NOT registered here: ch4_control_affine defers it straight to
% ch3_control_affine (M.m/V.m/G.m), so there is exactly one code path that
% ever calls the nominal dynamics, matching Chapter 3 bit for bit.
%
% Registered cases (tag = scale<SSS>, SSS = round(100*s)):
%   0.5  ->  scale050   ch4_test_model's ||Delta2||<1 feasibility-limit case
%   0.7  ->  scale070   Case III
%   1.5  ->  scale150   Case II
%   3.0  ->  scale300   Case IV
%
% To cover a new scale: run rabbit_generate_case_dynamics(s) once (writes
% M_scaleSSS.m / V_scaleSSS.m / G_scaleSSS.m into Dynamics/), then add s to
% CASES below. This function refuses an unregistered scale rather than
% silently falling back to s*M(q) -- that fallback is exactly the shortcut
% this design replaced.
%
% Inputs
%   s      : mass_scale; must match a registered case to within 1e-9
%   q      : nq x 1 joint position
%   dq     : nq x 1 joint velocity; omit (or pass []) if Vv is not needed
%
% Outputs
%   Mq, Vv, Gv : the case's mass matrix, Coriolis vector, gravity vector
%   tag        : the matched case tag, for diagnostics
%
% See also RABBIT_GENERATE_CASE_DYNAMICS, CH4_CONTROL_AFFINE, CH4_IMPACT.

CASES = [0.5 0.7 1.5 3.0];

[err, idx] = min(abs(CASES - s)); %#ok<ASGLU>
if abs(CASES(idx) - s) > 1e-9
    error('ch4_case_dynamics:unregistered', ...
          ['mass_scale %.6g has no pre-generated dynamics. Registered ' ...
           'cases: %s.\nRun rabbit_generate_case_dynamics(%.6g) once, then ' ...
           'add %.6g to CASES in ch4_case_dynamics.'], ...
          s, mat2str(CASES), s, s);
end

tag  = sprintf('scale%03.0f', round(CASES(idx)*100));
Mfun = str2func(['M_' tag]);
Gfun = str2func(['G_' tag]);

Mq = Mfun(q);
Gv = Gfun(q);

if nargin > 2 && ~isempty(dq)
    Vfun = str2func(['V_' tag]);
    Vv   = Vfun([q; dq]);
else
    Vv = [];
end

end
