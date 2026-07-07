function x_reset = rabbit_reset_map(x_post_impact)
    % 1. Relabel legs
    x_relabeled = relabel_state(x_post_impact);
    
    % 2. Calculate the position of the NEW stance foot in the world frame
    % We use the relabeled q (first 7 elements)
    q_rel = x_relabeled(1:7);
    P_new_stance = P_st(q_rel); 
    
    % 3. Apply the shift to the hip coordinates
    % q(1) = x_hip, q(2) = y_hip
    x_reset = x_relabeled;
%     x_reset(1) = x_reset(1) - P_new_stance(1);
    x_reset(2) = x_reset(2) - P_new_stance(2);
end
