function L = ch4_link_params()
%CH4_LINK_PARAMS  Per-body inertial parameters of today's 74 kg model.
%
%   L = ch4_link_params()
%
% The five bodies in the order ch4_link_dynamics uses them:
%
%   1 torso   2 stance thigh   3 stance shank   4 swing thigh   5 swing shank
%
% Fields (1 x 5 struct array):
%   .name
%   .m   mass [kg]
%   .a   COM offset ACROSS the link, along its local x axis [m]
%   .b   COM offset ALONG the link, along its local y axis [m] (the torso's
%        points up from the hip, so its b is negative in its own frame)
%   .J   planar moment of inertia about the COM, Izz [kg m^2]
%   .l   length from the proximal joint to the next joint or the foot [m];
%        geometry, not inertia, and never perturbed
%
% These are exactly Mass_Properties() of rabbit_energy_model_generalized_
% Lagrange.m as regenerated on 2026-09-02 (torso 47 kg, thighs 10 kg, shanks
% 3.5 kg; Izz 0.94 / 0.2 / 0.07), and ch4_test_model checks that
% ch4_link_dynamics(q, dq, ch4_link_params()) reproduces M.m / V.m / G.m.
% Ixx and Iyy of that file do not enter a planar model: rotation is about the
% out-of-plane axis only.
%
% See also CH4_LINK_DYNAMICS, CH4_UNCERTAINTY_SET.

L = struct('name', {'torso', 'thigh_st', 'shank_st', 'thigh_sw', 'shank_sw'}, ...
           'm',    {47,      10,         3.5,        10,         3.5}, ...
           'a',    {0,       0,          0,          0,          0}, ...
           'b',    {-0.375,  0.25,       0.25,       0.25,       0.25}, ...
           'J',    {0.94,    0.2,        0.07,       0.2,        0.07}, ...
           'l',    {0.75,    0.5,        0.5,        0.5,        0.5});

end
