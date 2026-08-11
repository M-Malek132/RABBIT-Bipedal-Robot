function stone = ch6_resolve_stone(sp, l_min, l_max)
%CH6_RESOLVE_STONE  Turn a terrain description into the barrier's live foothold.
%
%   stone = ch6_resolve_stone(sp, l_min, l_max)
%
% p.stones (sp here) describes the TERRAIN -- how stones move, how R1 and R2 are
% chosen. p.stone describes the ONE foothold the barrier is looking at during
% the current step, with every radius already a number. This is the conversion,
% and it happens once per step rather than once per QP solve.
%
% ----------------------------------------------------- why R2 is frozen per step
% With R2_mode = 'clearance', R2 is derived from l_min. Sections 6.2's stones
% move, so l_min changes continuously -- and if R2 followed it, circle O2 would
% change shape as well as position, the barrier's time derivatives would pick up
% two more terms, and the swing-foot clearance the robot is being asked for
% would drift mid-step for no physical reason. R2 encodes a DESIGN CHOICE about
% the swing trajectory. It is resolved here, from the foothold as it stands at
% the start of the step, and then held.
%
% l_min and l_max themselves are NOT frozen: they are stored as l_min0 / l_max0
% and ch6_stone_level adds the motion on top.
%
% Inputs
%   sp    : p.stones -- .R1 .R2_mode .R2_clearance .R2_fixed .motion
%           .v_stone .amp .freq
%   l_min : lower edge of this step's foothold window, at t = 0 [m]
%   l_max : upper edge [m]
%
% Output
%   stone : struct .l_min0 .l_max0 .R1 .R2 .clearance .motion .v_stone
%                  .amp .freq
%
% See also CH6_BAR_STONES, CH6_STONE_LEVEL, CH6_R2_FROM_CLEARANCE, CH6_SIMULATE.

if l_max < l_min
    error('ch6_resolve_stone:window', ...
          'Empty foothold window: l_max = %.4f < l_min = %.4f.', l_max, l_min);
end

switch lower(sp.R2_mode)
    case 'clearance'
        [R2, c_used] = ch6_R2_from_clearance(l_min, sp.R2_clearance);
    case 'fixed'
        R2     = sp.R2_fixed;
        c_used = sqrt(R2^2 + (l_min/2)^2) - R2;   % the clearance it implies
    otherwise
        error('ch6_resolve_stone:R2_mode', ...
              'Unknown R2_mode "%s" (expected clearance|fixed).', sp.R2_mode);
end

stone = struct('l_min0',    l_min, ...
               'l_max0',    l_max, ...
               'R1',        sp.R1, ...
               'R2',        R2, ...
               'clearance', c_used, ...
               'motion',    sp.motion, ...
               'v_stone',   sp.v_stone, ...
               'amp',       sp.amp, ...
               'freq',      sp.freq);

end
