function e = ch5_ecbf_gain(p, rb)
%CH5_ECBF_GAIN  Design the Exponential CBF gain Kb, eq (5.23)-(5.27), Thm 5.2.
%
%   e = ch5_ecbf_gain(p, rb)
%
% This is the function that makes Section 5.2 a DESIGN METHOD rather than an
% existence result. Applying VIOL to a relative-degree-rb constraint puts the
% barrier dynamics in controllable canonical form, (5.23):
%
%       etadot_b = Fb eta_b + Gb mu_b,     h = Cb eta_b
%
% -- a chain of rb integrators, and nothing about the plant survives into it.
% So the barrier is designed with POLE PLACEMENT, on a linear system, whatever
% the plant was.
%
% ----------------------------------------------------- from driving to bounding
% The step from stabilization to safety is (5.17) -> (5.18), and it is one
% character wide. To drive h to zero you would apply
%
%       mu_b = -Kb eta_b     =>   etadot_b =  Ab eta_b,   Ab = Fb - Gb Kb
%
% and get h(t) = Cb exp(Ab t) eta_b(0) -> 0. To keep h nonnegative you instead
% REQUIRE
%
%       mu_b >= -Kb eta_b    =>   etadot_b >= Ab eta_b               (5.25)
%
% and get h(t) >= Cb exp(Ab t) eta_b(0) >= 0. The equality that was a control
% law becomes an inequality that is a constraint, and the exponential that was
% a convergence rate becomes a LOWER ENVELOPE on the barrier. That is what the
% name Exponential CBF refers to.
%
% ------------------------------------------------ Theorem 5.2's two conditions
% Ab must be Hurwitz AND "TOTAL NEGATIVE" -- real negative eigenvalues, not
% merely negative real parts. Both are checked here, and the second is the one
% worth being strict about:
%
%   Theorem 5.1 is proved by applying Proposition 5.1 rb times down the chain
%   of sets C_rb -> C_rb-1 -> ... -> C_0, and each link uses the recursion
%   (5.30), y_i = ydot_(i-1) + p_i y_(i-1), with p_i a REAL positive number. A
%   complex conjugate pair is perfectly Hurwitz and gives an Ab with the right
%   characteristic polynomial, but there is no real p_i to build those
%   intermediate sets from, so the induction has no rungs and forward
%   invariance of C_0 does not follow. The exponential envelope would also
%   oscillate, which is precisely the un-smooth behaviour the chapter's
%   "smooth boundary for the safety set" is contrasting itself against.
%
% ------------------------------------------------------------- what is NOT here
% Corollary 5.2's condition p_i >= -ydot_(i-1)(x0)/y_(i-1)(x0) depends on the
% INITIAL CONDITION and so cannot be checked from p alone. It lives in
% ch5_ecbf_admissible, and ch5_simulate calls it before it starts integrating.
% This split is deliberate: the chapter is explicit that pole choice depending
% on x0 is a real limitation of the method, so the x0-dependent part is a
% separate, visible step rather than folded into the gain design.
%
% Inputs
%   p  : parameter struct (uses p.ecbf.poles)
%   rb : relative degree of the constraint
%
% Output
%   e : struct with
%         .poles .Kb (1 x rb) .Ab .Fb .Gb .Cb .rb
%         .eig_Ab .hurwitz .total_negative
%
% See also CH5_POLE_GAIN, CH5_ECBF_ADMISSIBLE, CH5_CTRL_ECBF_CLF_QP.

poles = p.ecbf.poles(:).';

if numel(poles) ~= rb
    error('ch5_ecbf_gain:poleCount', ...
          ['The constraint has relative degree %d, so p.ecbf.poles needs %d ' ...
           'entries (got %d: %s).'], rb, rb, numel(poles), mat2str(poles, 4));
end

Fb = diag(ones(rb-1,1), 1);            % (5.23)
Gb = [zeros(rb-1,1); 1];
Cb = [1, zeros(1, rb-1)];              % (5.11)

Kb = ch5_pole_gain(poles, rb);
Ab = Fb - Gb * Kb;

ev = eig(Ab);

% Tolerances are relative to the pole magnitudes: a pole set spanning 1.2 to
% 1.8 and one spanning 5 to 10 do not deserve the same absolute threshold.
scale = max(1, max(abs(poles)));

hurwitz = all(real(ev) < 0);
total_negative = hurwitz && all(abs(imag(ev)) < 1e-8 * scale);

if ~hurwitz
    error('ch5_ecbf_gain:notHurwitz', ...
          'Ab is not Hurwitz (max real eig %.3e) -- check p.ecbf.poles.', ...
          max(real(ev)));
end
if ~total_negative
    error('ch5_ecbf_gain:notTotalNegative', ...
          ['Ab is Hurwitz but NOT total negative (max |imag| eig %.3e). ' ...
           'Theorem 5.2 needs real negative eigenvalues: the proof of ' ...
           'Theorem 5.1 chains Proposition 5.1 through the real first-order ' ...
           'filters of (5.30), and a complex pair leaves that induction ' ...
           'without intermediate sets. Use real positive p.ecbf.poles.'], ...
          max(abs(imag(ev))));
end

e = struct('poles', poles, 'Kb', Kb, 'Ab', Ab, 'Fb', Fb, 'Gb', Gb, 'Cb', Cb, ...
           'rb', rb, 'eig_Ab', ev.', 'hurwitz', hurwitz, ...
           'total_negative', total_negative);

end
