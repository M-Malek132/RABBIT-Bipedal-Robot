function V = ch4_validity(sim, p)
%CH4_VALIDITY  How much of a Chapter-4 run a real robot could have walked.
%
%   V = ch4_validity(sim, p)
%
% The score is ch3_validity's, unchanged: a sample is invalid when the stance
% contact force lifts off (Fz <= 0) or asks for more than p.limits.mu_s of
% friction, and a run is valid up to its first invalid sample. What makes it a
% Chapter-4 score is the force it reads. ch4_step records .lambda from the TRUE
% model under the torque held there -- not the controller's nominal prediction,
% which is what its friction rows constrain (Remark 4.4) -- and, under a load
% redrawn every step, under the load that step carried.
%
% Inputs / output: as ch3_validity.
%
% See also CH3_VALIDITY, CH4_STEP, CH4_RUN_ENTRY.

V = ch3_validity(sim, p);

end
