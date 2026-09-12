function pc = ch4_run_params(p, controller, mass_scale, u_box)
%CH4_RUN_PARAMS  Parameters for ONE run of a Chapter-4 comparison.
%
%   pc = ch4_run_params(p, controller, mass_scale, u_box)
%   pc = ch4_run_params(p, controller, mass_scale, [])     keep p's boxes
%
% The single place that decides what a controller is told in a comparison run:
% which law, which model error it faces without knowing it, and how big its
% torque box is. ch4_compare_controllers builds the tables from it and
% ch4_animate builds the GIFs from it, so the picture cannot drift from the
% table. It did once: the animation left limits.u_max at the Chapter-3 default
% of 120 Nm while the sweep ran its per-case boxes, so the GIF raced a
% different controller than the one the table scored.
%
% ONLY THE LAWS WHOSE FORMULATION CARRIES A BOX ARE TOLD ABOUT ONE. The
% baselines are not -- that asymmetry is the experiment.
%
% THE CONTACT ROWS ARE LEFT AS p HAS THEM -- on, by default -- and that is a
% measured choice. (4.13) is the Chapter-3 stage-8 QP with the robust CLF row:
% torque box, friction cone and normal-force floor, the last two written on the
% NOMINAL model (Remark 4.4). Cutting it down to the torque box alone reads like
% the cleaner version of "robust CLF-QP with torque saturation", and it is not
% physical. Once eta ~= 0 the robust law carries a term of fixed magnitude
% D1/(1-D2) -- about 574 rad/s^2 on posture_195 -- along -LgV', whose
% direction flips whenever LgV = 2 eta' Peps G passes through zero however
% small eta is, so sampled at 1 kHz it chatters. With the
% floor removed that chatter pulled the TRUE stance foot into the ground for
% 17 / 38 / 42% of the samples in Cases I-III (min Fz -1621 / -3220 / -3543 N),
% which the bilateral stance model integrates without complaint. With the rows
% on, true Fz stayed >= 23 N in every case and the robust law still walked
% all three steps where the baselines fell. 'l1_con' drops the rows itself
% (ch4_ctrl_l1): mu2 is applied outside its QP, so they could not bound the
% realized torque anyway.
%
% Inputs
%   p          : parameter struct
%   controller : p.controller name for this run
%   mass_scale : the TRUE model's mass scale for this run
%   u_box      : torque box [Nm] for the law that carries one, applied to both
%                p.limits.u_max (robust laws) and p.l1.u_max ('l1_con'); []
%                leaves both as p already has them
%
% Output
%   pc : parameter struct for ch4_simulate / ch4_forces
%
% See also CH4_COMPARE_CONTROLLERS, CH4_ANIMATE, CH4_CTRL_RCLF_QP, CH4_CTRL_L1.

pc = p;
pc.controller             = controller;
pc.uncertainty.mass_scale = mass_scale;

if ~isempty(u_box)
    pc.limits.u_max = u_box;
    pc.l1.u_max     = u_box;
end

pc.limits.enable.torque = any(strcmpi(controller, {'clfqp_con', 'rclfqp_con'}));

end
