function x_new = relabel_state(x)
    % RELABEL_STATE Swaps stance and swing leg indices.
    % x = [q; dq] where q = [x; y; q_hip; q_st_h; q_st_k; q_sw_h; q_sw_k]
    
    x_new = x;
    % Indices: 1:x, 2:y, 3:q_h, 4:q_st_h, 5:q_st_k, 6:q_sw_h, 7:q_sw_k
    % Swap stance (4,5) with swing (6,7)
    idx = [4, 5, 6, 7];
    x_new([idx(3), idx(4), idx(1), idx(2)]) = x(idx);
    
    % Swap velocities (assuming dq follows same index pattern)
    idx_dq = idx + 7;
    x_new([idx_dq(3), idx_dq(4), idx_dq(1), idx_dq(2)]) = x(idx_dq);
end
