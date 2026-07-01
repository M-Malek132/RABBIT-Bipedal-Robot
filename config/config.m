function varargout = config(varargin)
    persistent data;
    
    % Handle the empty call case
    if nargin == 0
        varargout{1} = data;
        return;
    end
    
    action = varargin{1};
    switch action
        case 'init'

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
            
            % 3. Controller Parameters
            p.ctrl = init_bspline_params();
            
            data = p;
            fprintf('Configuration initialized: System constants and stones loaded.\n');
            varargout{1} = data;
        case 'get'
            varargout{1} = data;
    end
end