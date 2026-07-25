function x_out = ch3_relabel(x_in, p)
%CH3_RELABEL  Swap stance/swing leg coordinates (the "S" map of the reset).
%
%   x_out = ch3_relabel(x_in, p)
%
% The foot that just landed IS the new stance foot, so the two legs exchange
% roles.  RABBIT's legs are geometrically identical, so the relabeling is
% exactly an index swap -- no physical remapping of angles is needed.
%
%   q  = [ px pz qt | q1 q2 | q3 q4 ]   idx 1..7   (stance | swing)
%   dq = [ same order ]                 idx 8..14
%
% Only the four leg joints (4,5 <-> 6,7) and their velocities (11,12 <-> 13,14)
% move.  The floating base (px, pz, qt) and its velocity are untouched.
%
% See also CH3_IMPACT.

nq  = p.nq;
st  = [4 5];        % stance leg joints
sw  = [6 7];        % swing  leg joints

x_out = x_in;

x_out(st) = x_in(sw);
x_out(sw) = x_in(st);

x_out(st + nq) = x_in(sw + nq);
x_out(sw + nq) = x_in(st + nq);

end
