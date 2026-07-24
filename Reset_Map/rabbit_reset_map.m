function x_reset = rabbit_reset_map(x_post_impact)
%RABBIT_RESET_MAP  Post-impact reset: relabel legs + re-plant new stance foot.
%
%   x_reset = rabbit_reset_map(x_post_impact)
%
% The second half of the hybrid transition, applied AFTER rabbit_impact_map:
%   1. relabel_state swaps the stance/swing leg roles (the foot that just
%      landed is the new stance foot);
%   2. the vertical base coordinate pz is shifted so the NEW stance foot sits
%      exactly on the ground (z = 0), correcting any small residual foot height.
%
% Callers always compose the two as
%       x_reset = rabbit_reset_map(rabbit_impact_map(x_minus)).
%
% The HORIZONTAL base coordinate px is deliberately NOT shifted (line below is
% left commented): px is the translation gauge and must advance with the robot,
% so it accumulates across steps. Every periodicity check excludes px (index 1)
% for exactly this reason (see hzd_constraints / poincare_stability).
%
% Input  : x_post_impact 14x1 post-impact, pre-relabel state
% Output : x_reset       14x1 start-of-next-step state (legs relabeled,
%                        new stance foot on the ground)

    % 1. Relabel legs (swing foot just landed -> it is the new stance foot).
    x_relabeled = relabel_state(x_post_impact);
    
    % 2. Calculate the position of the NEW stance foot in the world frame
    % We use the relabeled q (first 7 elements)
    q_rel = x_relabeled(1:7);
    P_new_stance = P_st(q_rel); 
    
    % 3. Apply the shift to the hip coordinates
    % q(1) = x_hip, q(2) = y_hip
    x_reset = x_relabeled;
%     x_reset(1) = x_reset(1) - P_new_stance(1);

    % SIGN: foot_z = -pz - C, so driving foot_z to 0 needs pz_new = pz + foot_z.
    % This was a MINUS, which doubles the foot-height error instead of removing
    % it. Mostly masked in normal operation (the impact event fires when the
    % swing foot is already at z~0, so the correction is ~0 either way), but it
    % is wrong and disagreed with the identical projection in
    % make_initial_state.m and make_random_initial_state.m, which both use +.
    x_reset(2) = x_reset(2) + P_new_stance(2);
end
