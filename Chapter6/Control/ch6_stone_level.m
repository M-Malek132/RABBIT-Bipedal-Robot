function s = ch6_stone_level(stone, t)
%CH6_STONE_LEVEL  The foothold window [l_min(t), l_max(t)] and its derivatives.
%
%   s = ch6_stone_level(stone, t)
%
% Section 6.1.3 has a stone that sits still, so l_min and l_max are constants
% and everything below is zero. Section 6.2 has a stone that MOVES, and the
% barrier then needs the first and second time derivatives of the window as
% well as its value, because g depends on t explicitly and
%
%       gdot  = dg/dt + ...          gddot = d2g/dt2 + 2 d2g/dtdr rdot + ...
%
% Returning the derivatives from the same place that defines the motion is what
% keeps them consistent: a motion law added here without its derivatives fails
% loudly (the switch has no default branch) instead of quietly contributing a
% zero feed-forward.
%
% ------------------------------------------------------------------- freezing
% "These stepping stones move with time, stopping once a foot is placed on it."
% The stone therefore stops at IMPACT, not at some absolute time -- which is why
% t here is time SINCE THE START OF THE CURRENT STEP and why ch6_simulate
% re-seeds stone.l_min0 / l_max0 from the frozen values when it moves on. There
% is nothing to do at this level: the step ends, and the next step's stone
% starts its own clock.
%
% -------------------------------------------------------------------- motions
%   'static'      l(t) = l0
%   'linear'      l(t) = l0 + v t                        both edges drift
%   'sinusoidal'  l(t) = l0 + A sin(2 pi f t)            both edges oscillate
%
% Both edges move together, so the stone TRANSLATES and its size l_max - l_min
% is constant. A stone that changed size while moving would confound "can the
% controller track a moving target" with "can it hit a shrinking one", and
% Fig. 6.10 is about the first.
%
% Inputs
%   stone : struct with .l_min0 .l_max0 .motion and, as the motion needs them,
%           .v_stone .amp .freq
%   t     : time since the start of the current step [s]
%
% Output
%   s : struct .l_min .l_max .dl_min .dl_max .ddl_min .ddl_max
%
% See also CH6_BAR_STONES, CH6_TERRAIN, CH6_SIMULATE.

switch lower(stone.motion)

    case 'static'
        d  = 0;
        dd = 0;
        o  = 0;

    case 'linear'
        o  = stone.v_stone * t;
        d  = stone.v_stone;
        dd = 0;

    case 'sinusoidal'
        w  = 2*pi*stone.freq;
        o  =  stone.amp * sin(w*t);
        d  =  stone.amp * w * cos(w*t);
        dd = -stone.amp * w^2 * sin(w*t);

    otherwise
        error('ch6_stone_level:motion', ...
              'Unknown stone motion "%s" (expected static|linear|sinusoidal).', ...
              stone.motion);
end

s = struct('l_min',   stone.l_min0 + o, 'l_max',   stone.l_max0 + o, ...
           'dl_min',  d,                'dl_max',  d, ...
           'ddl_min', dd,               'ddl_max', dd);

end
