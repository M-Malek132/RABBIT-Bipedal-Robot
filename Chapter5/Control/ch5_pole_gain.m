function K = ch5_pole_gain(poles, r)
%CH5_POLE_GAIN  Pole-placement gain for a chain of r integrators.
%
%   K = ch5_pole_gain(poles, r)
%
% For the controllable canonical pair (5.23)
%
%       Fb = shift_r (ones on the superdiagonal),   Gb = e_r
%
% the closed loop Ab = Fb - Gb K has last row -K, and a companion matrix with
% last row [-k1 ... -kr] has characteristic polynomial
%
%       s^r + kr s^(r-1) + ... + k2 s + k1.
%
% Matching that to the desired  prod_i (s + p_i)  = s^r + a1 s^(r-1) + ... + ar
% gives the gains IN REVERSE ORDER:
%
%       K = [ar, a_{r-1}, ..., a1].
%
% The reversal is the whole content of this function and is easy to get
% backwards -- it produces a Hurwitz matrix either way, just with the wrong
% poles, so nothing downstream would complain. ch5_test_ecbf asserts
% eig(Fb - Gb K) == -poles rather than merely asserting stability.
%
% Inputs
%   poles : r positive reals; the closed-loop eigenvalues are their NEGATIVES
%   r     : chain length
%
% Output
%   K : 1 x r
%
% See also CH5_ECBF_GAIN, CH5_RES_CLF.

poles = poles(:).';

if numel(poles) ~= r
    error('ch5_pole_gain:count', ...
          'Need exactly %d poles for a relative-degree-%d chain (got %d).', ...
          r, r, numel(poles));
end
if any(~isreal(poles)) || any(poles <= 0)
    error('ch5_pole_gain:sign', ...
          ['Poles must be real and strictly positive; they are used as ' ...
           '-p_i. Got %s.'], mat2str(poles, 4));
end

a = poly(-poles);          % [1, a1, ..., ar]
K = fliplr(a(2:end));

end
