function [B, geom] = ch6_bar_obstacle(k, ob)
%CH6_BAR_OBSTACLE  The overhead-obstacle barriers of Section 6.1.2.
%
%   [B, geom] = ch6_bar_obstacle(k, ob)
%
% Two constraints on the head (torso top), in coordinates (l_H, h_H) measured
% from the stance foot, selected by ob.type:
%
%   'ceiling'  g_C = h_r - h_H(q) >= 0                                  (6.2)
%   'circle'   g_O = |(l_H,h_H) - (l_m, h_m + R1o)| - R1o >= 0          (6.5)
%
% ---------------------------------------------------- what the second one buys
% (6.2) costs head height for the WHOLE STEP whether or not the obstacle is
% anywhere near. (6.5) is the same requirement localised: the keep-out disc has
% its centre O1 one radius ABOVE the obstacle M at (l_m, h_m), so the disc rests
% on the obstacle, and staying outside it is a constraint that relaxes to
% nothing as the head moves away in l_H. Fig. 6.1's green circle is that disc.
%
% Both are drawn in Fig. 6.4 and the point of showing them together is that the
% second lets the robot walk taller everywhere except under the obstacle -- the
% same trade the stepping-stone circles make against a naive box on l_f.
%
% ------------------------------------------------------- a note on R1o's role
% R1o is NOT the size of the obstacle. The obstacle is the point M; R1o is how
% far ahead, in head travel, the constraint starts to bite. A small R1o gives a
% constraint that switches on abruptly as the head arrives, which the barrier
% must then answer with a large gddot; a large R1o approaches the flat ceiling.
% It is the same conservatism dial as R1 in Fig. 6.3 and p.width.l3 in Fig 6.12,
% and in all three cases the chapter leaves it to be picked by trial.
%
% Inputs
%   k  : head kinematics from ch6_kin(x, aux, p, 'head')
%   ob : p.obstacle struct -- .type .h_r .l_m .h_m .R1o
%
% Outputs
%   B    : 1x1 struct of barrier terms from ch6_bar_lift
%   geom : struct .l_H .h_H and the obstacle description, for plots
%
% See also CH6_KIN, CH6_BAR_LIFT, CH6_BAR_CIRCLE, CH6_BAR_AFFINE.

switch lower(ob.type)

    case 'ceiling'
        % g = h_r - h_H = h_r + [0 -1] r
        G     = ch6_bar_affine(k.r, [0; -1], ob.h_r);
        label = 'hH <= hr';

    case 'circle'
        cO    = [ob.l_m; ob.h_m + ob.R1o];
        G     = ch6_bar_circle(k.r, cO, ob.R1o, 'outside');
        label = 'head outside O1';

    otherwise
        error('ch6_bar_obstacle:type', ...
              'Unknown obstacle type "%s" (expected ceiling|circle).', ob.type);
end

B = ch6_bar_lift(G, k, label);

geom = struct('l_H', k.r(1), 'h_H', k.r(2), 'type', ob.type, ...
              'h_r', ob.h_r, 'l_m', ob.l_m, 'h_m', ob.h_m, 'R1o', ob.R1o);

end
