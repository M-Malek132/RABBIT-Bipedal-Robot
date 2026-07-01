function val = config(request)
    % CONFIG: Acts as a central database for all system parameters.
    % Uses existing parameters() and init_bspline_params() functions.
    
    persistent data;
    
    if nargin < 1, request = 'get'; end
    
    switch request
        case 'init'
            % 1. Get physical params from your existing function
            p = parameters(); 
            
            % 2. Generate Stones (Persistent across calls)
            num_stones = 20;
            stones = zeros(num_stones + 2, 2);
            stones(1, :) = [-0.5, 0.5];
            current_x = 0.5;
            for i = 1:num_stones
                gap = 0.05 + 0.15 * rand();
                width = 0.3 + 0.3 * rand();
                stones(i + 1, :) = [current_x + gap, current_x + gap + width];
                current_x = stones(i + 1, 2);
            end
            stones(end, :) = [current_x + 0.1, current_x + 6.0];
            p.stones = stones;
            
            % 3. Embed controller params from your existing function
            p.ctrl = init_bspline_params();
            
            data = p;
            fprintf('Configuration initialized: Robot parameters, B-splines, and Stones loaded.\n');
            
        case 'get'
            val = data;
            
        case 'stones'
            val = data.stones;
    end
end
