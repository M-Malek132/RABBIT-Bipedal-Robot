function b = ch5_barrier(x, p)
%CH5_BARRIER  The constraint function h and everything both CBF forms need.
%
%   b = ch5_barrier(x, p)
%
% ONE SIGN CONVENTION, EVERYWHERE IN CHAPTER 5:
%
%       SAFE SET     C = { x : h(x) >= 0 }                            (5.1)
%
% An upper bound q <= q_max becomes h = q_max - q; a lower bound q >= q_min
% becomes h = q - q_min. Nothing downstream needs to know which it was.
%
% ------------------------------------------------- the two barriers, separated
% The thesis writes "B(x)" for two different objects and this is where they get
% told apart, because they have OPPOSITE senses and no dimension check would
% catch a swap.
%
%   SECTION 5.1, RECIPROCAL (5.6):   B_rec = 1/h,  blows up at the boundary,
%     condition (5.3)  Bdot <= gamma/B.  Substituting B = 1/h and B_recdot =
%     -hdot/h^2 turns that into the far more readable
%
%           hdot >= -gamma h^3                                        (*)
%
%     which is returned as .rec_row so callers never re-derive it. Note (*) is
%     a condition on hdot, so IT ONLY CONSTRAINS u WHEN h HAS RELATIVE DEGREE 1.
%
%   SECTION 5.2, EXPONENTIAL (Def 5.1): h itself is the barrier, and the
%     condition is placed on its rb-th derivative,
%
%           h^(rb) >= -Kb eta_b,   eta_b = [h; hdot; ...; h^(rb-1)]   (5.10)
%
%     which constrains u at ANY relative degree, because h^(rb) is by
%     definition the first derivative in which u appears.
%
% ------------------------------------------------- the failure, made returnable
% .Lgh is L_g h. For a relative-degree-rb constraint with rb > 1 it is
% IDENTICALLY ZERO -- that is what relative degree > 1 means. So the Section
% 5.1 barrier row degenerates to
%
%           0 * u  >=  -gamma h^3 - L_f h
%
% a statement with no control in it, which the QP can neither satisfy nor
% violate. This is not a numerical fragility to be guarded against; it is the
% chapter's motivating limitation, and .rec_controllable reports it as data so
% ch5_ctrl_cbf_clf_qp can fail loudly and ch5_test_qp can assert on it.
%
% -------------------------------------------------------------- constraints
%   springmass  'x3_max'  x3 <= level        relative degree 6
%               'v1_max'  xdot1 <= level     relative degree 1
%   pendulum    'py_min'  py2 >= level       relative degree 4
%
% Output
%   b : struct with
%         .type .level .rb .q .sense .label
%         .h              scalar, >= 0 is safe
%         .eta_b          rb x 1,  [h; hdot; ...; h^(rb-1)]      (5.10)
%         .Lfrb           scalar,  L_f^rb h
%         .LgLfrb1        1 x nu,  L_g L_f^(rb-1) h
%         .Lfh .Lgh       the relative-degree-1 pieces, for Section 5.1
%         .B_rec .LfB_rec .LgB_rec    the (5.6) reciprocal barrier
%         .rec_controllable   false when Lgh == 0, i.e. rb > 1
%         .rec_defined        false when h <= 0 (1/h is not a barrier there)
%
% See also CH5_LIE_ROWS, CH5_ECBF_GAIN, CH5_CTRL_CBF_CLF_QP, CH5_CTRL_ECBF_CLF_QP.

x  = x(:);
cs = p.constraint;
nu = p.sys.nu;

switch lower(cs.type)

    %% ------------------------------------------------ x3 <= level  (rb = 6)
    case 'x3_max'
        sm = ch5_springmass(p);
        b  = linear_barrier(x, sm.A, sm.B, -[0 0 1 0 0 0], cs.value, sm.n);
        b.q     = x(3);
        b.sense = 'upper';
        b.label = 'x_3';

    %% --------------------------------------------- xdot1 <= level  (rb = 1)
    case 'v1_max'
        sm = ch5_springmass(p);
        b  = linear_barrier(x, sm.A, sm.B, -[0 0 0 1 0 0], cs.value, sm.n);
        b.q     = x(4);
        b.sense = 'upper';
        b.label = '$\dot x_1$';

    %% -------------------------------------------- py2 >= level  (rb = 4)
    case 'py_min'
        pv = p.plant.pv;
        [py, py1, py2, py3, py4f, py4g] = ch5_pend_py_lie(x, pv);

        b = struct();
        b.rb       = 4;
        b.h        = py - cs.value;
        b.eta_b    = [b.h; py1; py2; py3];
        b.Lfrb     = py4f;
        b.LgLfrb1  = py4g(:).';
        b.Lfh      = py1;
        b.Lgh      = zeros(1, nu);     % relative degree 4, so identically zero
        b.q        = py;
        b.sense    = 'lower';
        b.label    = 'p^y_2';

    otherwise
        error('ch5_barrier:unknownConstraint', ...
              'Unknown constraint type "%s".', cs.type);
end

b.type  = cs.type;
b.level = cs.value;

%% ------------------------------------- the Section 5.1 reciprocal barrier
% Built unconditionally so its degeneracy is visible rather than avoided.
b.rec_controllable = norm(b.Lgh, inf) > 1e-12;
b.rec_defined      = b.h > 0;

if b.rec_defined
    b.B_rec    = 1 / b.h;
    b.LfB_rec  = -b.Lfh / b.h^2;
    b.LgB_rec  = -b.Lgh / b.h^2;
else
    % 1/h at h <= 0 is not a barrier -- it is negative or infinite, and (5.2)'s
    % sandwich between class-K functions of h has already failed. Report rather
    % than return a number that looks usable.
    b.B_rec   = inf;
    b.LfB_rec = nan;
    b.LgB_rec = nan(1, nu);
end

end

% ---------------------------------------------------------------------------
function b = linear_barrier(x, A, B, c, level, n)
%LINEAR_BARRIER  h(x) = c x + level for a linear plant, exact at every order.
%
% c is the row such that h = c x + level. For an upper bound q = e_i' x <= level
% that is c = -e_i', which puts the level in as a plain +level and keeps the
% "h >= 0 is safe" convention without a second sign to track.

L  = ch5_lie_rows(A, B, c, n);
rb = L.rd;

if isnan(rb)
    error('ch5_barrier:noRelativeDegree', ...
          ['The control never reaches this constraint within %d ' ...
           'derivatives -- it is uncontrollable, not high relative degree.'], n);
end

b = struct();
b.rb      = rb;
b.eta_b   = L.rows(1:rb, :) * x;
b.eta_b(1) = b.eta_b(1) + level;         % only h itself carries the offset
b.h       = b.eta_b(1);
b.Lfrb    = L.rows(rb+1, :) * x;
b.LgLfrb1 = L.gcoef(rb, :);              % c A^(rb-1) B
b.Lfh     = L.rows(2, :) * x;
b.Lgh     = L.gcoef(1, :);               % c B, zero unless rb == 1

end
