function T = ch6_terrain(n, range, stone_sz, seed)
%CH6_TERRAIN  A randomly generated set of discrete footholds (Fig. 6.2, 6.5).
%
%   T = ch6_terrain(n, range, stone_sz)
%   T = ch6_terrain(n, range, stone_sz, seed)
%
% n stones, each a window [l_min, l_max] of width stone_sz whose CENTRE l_d is
% drawn uniformly from `range`. That is the thesis's setup for both Fig. 6.5
% ("desired step lengths ... chosen randomly in the range 0.25 m to 0.6 m") and
% Table 6.1 ("10 randomly placed stepping stones with a stone size of 5 cm").
%
% The window is centred on l_d and NOT clipped to `range`: a stone whose centre
% is at the top of the range genuinely extends past it, and clipping would
% quietly shrink the hardest stones into easier ones -- which is the one
% distortion that would flatter the controller in exactly the cells of Table 6.1
% where it is under most pressure.
%
% ---------------------------------------------------------------- reproducibility
% Pass `seed` and the terrain is deterministic. Table 6.1 is a Monte-Carlo
% estimate over 100-ish trials whose whole content is a percentage, so a run
% that cannot be reproduced is a number nobody can check. ch6_table61 seeds
% every trial from p.mc.seed and the trial index, so an individual failing trial
% can be re-run on its own.
%
% Inputs
%   n        : number of stones
%   range    : [l_lo l_hi], bounds on the desired step length l_d [m]
%   stone_sz : l_max - l_min [m]
%   seed     : optional rng seed
%
% Output
%   T : 1 x n struct array with .l_d .l_min .l_max
%
% See also CH6_SIMULATE, CH6_TABLE61, CH6_RESOLVE_STONE.

if nargin >= 4 && ~isempty(seed)
    rs = RandStream('twister', 'Seed', seed);
    u  = rand(rs, 1, n);
else
    u  = rand(1, n);
end

l_d = range(1) + (range(2) - range(1)) * u;
h   = stone_sz / 2;

T = struct('l_d',   num2cell(l_d), ...
           'l_min', num2cell(l_d - h), ...
           'l_max', num2cell(l_d + h));

end
