function params = merge_params()
% Merge robot physical parameters with controller parameters
% Returns a single params structure for the simulator

    % Get robot physical parameters
    params = parameters();
    
    % Get controller parameters
    ctrl = init_bspline_params();
    
    % Embed controller params into params structure
    params.ctrl = ctrl;
    
    fprintf('Parameters merged:\n');
    fprintf('  Robot: mT=%.1f, m1=%.1f, m2=%.1f, l1=%.2f, l2=%.2f\n', ...
        params.mT, params.m1, params.m2, params.l1, params.l2);
    fprintf('  B-spline: n=%d, p=%d, CP=%d\n', ...
        ctrl.n, ctrl.p, ctrl.n+1);
    fprintf('  Phase: [%.2f, %.2f]\n', ctrl.theta0, ctrl.thetaf);
    
end