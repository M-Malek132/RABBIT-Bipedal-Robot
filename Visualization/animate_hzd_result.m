function animate_hzd_result(z_opt, p, nSteps)
    if nargin < 3, nSteps = 10; end

    [coeffs, x_current] = unpack_z(z_opt, p);

    t_all = [];
    x_all = [];
    time_offset = 0;

    for step = 1:nSteps
        [t_step, x_step] = simulate_hzd_gait_full(coeffs, x_current, p);

        if isempty(t_all)
            t_all = t_step;
            x_all = x_step;
        else
            t_all = [t_all; t_step(2:end) + time_offset];
            x_all = [x_all; x_step(2:end, :)];
        end

        x_minus        = x_step(end, :)';
        x_after_impact = rabbit_impact_map(x_minus);
        x_current      = rabbit_reset_map(x_after_impact);

        time_offset = t_all(end);
    end

    animate_rabbit(x_all);
end