function sig = ch3_model_signature(p)
%CH3_MODEL_SIGNATURE  A fingerprint of the dynamics a gait is solved on.
%
%   sig = ch3_model_signature()      today's model (M.m / V.m / G.m)
%   sig = ch3_model_signature(p)     the model p selects (p.model_blend)
%
% WHY. M.m, V.m and G.m are global generated files. A gait saved before they
% were regenerated still loads, unpacks and simulates -- as a reference that is
% no longer an orbit of the robot. That happened on 2026-09-02 (30 kg -> 74 kg,
% non-uniformly): five stored gaits kept loading, none was still an orbit, and
% the old reference gait was UNSTABLE on the new model (delta^2 0.904, rho
% 1.18). Nothing in the files could say so. ch3_col_solve now stores this
% signature with every solution (out.model_sig, and out.p.model_sig), and
% ch3_model_check compares it with the model on the path when a gait is loaded.
%
% WHAT IS FINGERPRINTED. The generated functions themselves, evaluated at three
% fixed states: M, V and G through ch3_mvg (so a blended model is told apart
% from today's), and the kinematics P_st, P_sw, J_st. Evaluating the dynamics
% rather than hashing a parameter file catches the failure that actually
% occurred -- files regenerated from edited parameters -- and also a hand edit
% or a partial regeneration. Hashing the .m text would not survive a harmless
% reformat or a different Symbolic Math Toolbox version; comparing numbers
% does.
%
% COMPARED BY VALUE, NOT BY DIGEST. Floating-point results can differ in the
% last bits between machines (this project runs on two), so ch3_model_check
% compares .vals to a relative tolerance. .digest is only a short label for
% logs and reports, built from the values rounded to 8 significant digits.
%
% Output
%   sig : struct
%     .vals     column vector of the evaluated quantities
%     .digest   12-hex-digit label
%     .mass     total mass implied by gravity, -G(2)/g0 [kg]
%     .blend    p.model_blend ([] = today's model)
%     .version  layout version of .vals (bump if the fingerprint changes)
%     .created  when it was computed
%
% See also CH3_MODEL_CHECK, CH3_MVG, CH3_COL_SOLVE, CH3_STAMP_GAITS.

if nargin < 1 || isempty(p)
    p = struct('model_blend', [], 'g0', 9.8062);
end
if ~isfield(p, 'g0'), p.g0 = 9.8062; end
blend = [];
if isfield(p, 'model_blend'), blend = p.model_blend; end

% Three fixed, generic states -- not random, so the fingerprint does not
% depend on the random-number state of whoever calls it.
Q  = [ 0.10 -0.95  0.05 -0.30  0.40  0.20  0.30; ...
      -0.20 -0.90 -0.10  0.25  0.10 -0.45  0.70; ...
       0.35 -0.85  0.20 -0.60  0.90  0.50  0.05 ].';
DQ = [ 0.50 -0.20  0.30 -1.10  0.70  0.90 -0.40; ...
      -1.00  0.40 -0.60  0.80 -1.50  0.30  1.20; ...
       0.20  0.10  0.90 -0.30  0.60 -1.40  0.50 ].';

vals = [];
for k = 1:size(Q, 2)
    q = Q(:, k);  dq = DQ(:, k);
    [Mq, Vv, Gv] = ch3_mvg(q, dq, p);
    vals = [vals; Mq(:); Vv(:); Gv(:); P_st(q); P_sw(q); ...
            reshape(J_st(q), [], 1)];               %#ok<AGROW>
end

[~, ~, G0] = ch3_mvg(Q(:, 1), [], p);

sig = struct('vals', vals, 'digest', digest(vals), ...
             'mass', -G0(2) / p.g0, 'blend', blend, 'version', 1, ...
             'created', datestr(now, 'yyyy-mm-dd HH:MM:SS')); %#ok<TNOW1,DATST>

end

% ---------------------------------------------------------------------------
function h = digest(v)
%DIGEST  Short hex label of the values rounded to 8 significant digits.
txt = sprintf('%.7e,', round(v, 8, 'significant'));
try
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(int8(double(txt)));
    b  = typecast(md.digest(), 'uint8');
    h  = sprintf('%02x', b(1:6));
catch
    % No JVM (matlab -nojvm): a polynomial checksum of the same text, in
    % doubles (exact: acc*131 + 127 stays far below 2^53). Not MD5, so the
    % label differs from a JVM session's -- which is harmless, since the check
    % compares .vals, never the label.
    acc = 0;
    for c = double(txt)
        acc = mod(acc * 131 + c, 4294967291);
    end
    h = sprintf('%08x', acc);
end
end
