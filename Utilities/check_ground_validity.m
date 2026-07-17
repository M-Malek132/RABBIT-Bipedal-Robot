function report = check_ground_validity(q, dq, tol)
    if nargin < 3, tol = 1e-3; end
    pts = get_body_points(q);

    % World-frame Z (from get_body_points/Tt/P_st) is UP-positive:
    % ground = 0, above ground = z > 0. A point with z < -tol has gone
    % BELOW ground (penetrating). max_penetration is the worst (largest)
    % depth below ground across the sampled points; negative means
    % everything is comfortably above ground by that margin.
    fields = {'torso_top','stance_knee','swing_knee','swing_foot'};
    report.max_penetration = -inf;
    report.violations = {};
    for i = 1:numel(fields)
        h = pts.(fields{i})(2);
        if h < -tol
            report.violations{end+1} = sprintf('%s penetrating ground: z=%.4f', fields{i}, h);
        end
        report.max_penetration = max(report.max_penetration, -h);
    end

    report.stance_foot_height = pts.stance_foot(2);
    if abs(report.stance_foot_height) > tol
        report.violations{end+1} = sprintf('stance foot NOT on ground: z=%.4f', report.stance_foot_height);
    end

    if nargin >= 2 && ~isempty(dq)
        J = J_st(q);
        v = J*dq;
        report.stance_foot_velocity = norm(v);
        if norm(v) > 1e-2
            report.violations{end+1} = sprintf('stance foot slipping: |v|=%.4f', norm(v));
        end
    end

    report.is_valid = isempty(report.violations);
end