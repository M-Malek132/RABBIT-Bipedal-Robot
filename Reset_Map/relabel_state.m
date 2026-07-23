function x_new = relabel_state(x)
%RELABEL_STATE  Swap stance and swing leg coordinates after an impact.
%
%   x_new = relabel_state(x)
%
% Immediately after the swing foot lands it becomes the new STANCE foot, so the
% roles of the two legs exchange. Because RABBIT's legs are geometrically
% identical, this relabeling is exactly an index swap of the leg joints (and
% their velocities); no physical remapping of angles is needed. This is the
% "S : x -> x" coordinate swap of Westervelt et al. (Sec. 3.4.2), applied here
% on the POST-impact state (see rabbit_reset_map, which calls this then
% re-plants the new stance foot on the ground).
%
% State layout (14x1), x = [q; dq]:
%   q  = [ px  pz  qt  q_st_hip  q_st_knee  q_sw_hip  q_sw_knee ]   idx 1..7
%   dq = [ ...same order... ]                                        idx 8..14
% The floating-base coordinates (px, pz, qt) and their velocities are NOT
% touched; only the four leg joints (idx 4,5 <-> 6,7) and the matching
% velocities (idx 11,12 <-> 13,14) are swapped.
%
% Input  : x     14x1 post-impact state (pre-relabel)
% Output : x_new 14x1 same state with stance/swing leg roles exchanged

    x_new = x;
    % Swap the leg joint POSITIONS: stance (4,5) <-> swing (6,7).
    idx = [4, 5, 6, 7];
    x_new([idx(3), idx(4), idx(1), idx(2)]) = x(idx);

    % Swap the matching leg joint VELOCITIES (same pattern, offset by nq = 7).
    idx_dq = idx + 7;
    x_new([idx_dq(3), idx_dq(4), idx_dq(1), idx_dq(2)]) = x(idx_dq);
end
