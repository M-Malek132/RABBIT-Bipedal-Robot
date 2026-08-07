function [B, geom, k] = ch6_barrier(x, aux, p, t)
%CH6_BARRIER  Single dispatch point for the Chapter-6 position constraints.
%
%   [B, geom, k] = ch6_barrier(x, aux, p, t)
%
% p.cbf.problem picks which of the chapter's constraint sets is active:
%
%   'stones'    Sections 6.1.3 / 6.2 / 6.4 -- footstep placement, 2 barriers
%   'obstacle'  Section 6.1.2 -- overhead obstacle, 1 barrier
%   'none'      no barrier rows; the QP degenerates to Chapter 3 stage 8
%
% Everything downstream -- the QP, the simulation, the report, the plots --
% takes a struct ARRAY of barriers and never asks which problem produced them,
% so adding a constraint set is adding a case here and nothing else. That is
% Remark 6.1 taken literally: "with this methodology, and different design for
% barrier constraints, we can further apply this approach for a variety of
% additional constraints".
%
% 'none' exists because Table 6.1 needs it. Controller I in (6.27) is the gait
% library with no CBF, and running it through the same code path as the others
% -- same integrator, same sampling, same torque box -- is what makes the
% columns of that table comparable.
%
% Inputs
%   x   : 14x1 state
%   aux : third output of ch3_control_affine, or [] to compute it
%   p   : parameter struct; for 'stones' it must carry p.stone, the resolved
%         per-step foothold (ch6_simulate sets it)
%   t   : time since the start of the current step [s]
%
% Outputs
%   B    : 1 x nb struct array from ch6_bar_lift (possibly empty)
%   geom : problem-specific geometry for plots and reports
%   k    : the kinematics struct, so callers reuse the one evaluation
%
% See also CH6_BAR_STONES, CH6_BAR_OBSTACLE, CH6_CBF_ROW, CH6_CONTROL.

if nargin < 4 || isempty(t), t = 0; end

switch lower(p.cbf.problem)

    case 'stones'
        k = ch6_kin(x, aux, p, 'swing');
        [B, geom] = ch6_bar_stones(k, p.stone, t);

    case 'obstacle'
        k = ch6_kin(x, aux, p, 'head');
        [B, geom] = ch6_bar_obstacle(k, p.obstacle);

    case 'none'
        k    = struct('r', [NaN; NaN]);
        B    = ch6_bar_lift(ch6_bar_affine([0;0], [0;0], 1), ...
                            struct('rdot', [0;0], 'Lf2r', [0;0], ...
                                   'LgLfr', zeros(2, p.nu)), 'none');
        B    = B([]);                 % empty 1x0 array of the right type
        geom = struct();

    otherwise
        error('ch6_barrier:problem', ...
              ['Unknown p.cbf.problem "%s" ' ...
               '(expected stones|obstacle|none).'], p.cbf.problem);
end

end
