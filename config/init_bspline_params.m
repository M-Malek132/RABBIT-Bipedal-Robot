function ctrl = init_bspline_params()
    ctrl.n = 7; ctrl.p = 3;
    ctrl.theta0 = 0.05; ctrl.thetaf = 0.35;
    
    qa_start = [-0.3; 0.6; -1.0; 0.6];
    qa_end   = [ 0.3; -0.3; 0.3; -0.3];
    
    ctrl.ControlPoints = zeros(ctrl.n+1, 4);
    for j = 1:4
        ctrl.ControlPoints(1, j) = qa_start(j);
        ctrl.ControlPoints(end, j) = qa_end(j);
        for k = 1:ctrl.n-1
            s_k = k / ctrl.n;
            ctrl.ControlPoints(k+1, j) = qa_start(j) + s_k * (qa_end(j) - qa_start(j));
        end
    end
    
    ctrl.Kp = diag([300, 300, 200, 200]);
    ctrl.Kd = diag([40,  40,  25,  25 ]);
end
