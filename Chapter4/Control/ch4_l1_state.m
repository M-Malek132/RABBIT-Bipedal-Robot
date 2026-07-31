function out = ch4_l1_state(action, p, varargin)
%CH4_L1_STATE  Layout, initialization and reset of the L1 controller state.
%
%   n   = ch4_l1_state('dim',    p)
%   xi  = ch4_l1_state('init',   p, eta0)
%   s   = ch4_l1_state('unpack', p, xi)
%   xi  = ch4_l1_state('pack',   p, s)
%   xi  = ch4_l1_state('reset',  p, xi, eta_plus)
%
% Unlike every controller in Chapter 3, the L1 controller HAS MEMORY. It is not
% a function of x; it is a dynamical system driven by x, and the simulation has
% to integrate it alongside the robot. This file defines that state once so the
% ODE right-hand side, the step integrator and the analysis code all agree on
% what the 5*ny entries mean.
%
% THE STATE, in order (5*ny = 20 entries for RABBIT):
%
%   eta_hat    2ny   state predictor (4.19). Its job is to reproduce the real
%                    transverse dynamics; the gap eta_tilde = eta_hat - eta is
%                    the ONLY signal the adaptation gets to see.
%   alpha_hat   ny   estimate of the state-proportional part of theta (4.17)
%   beta_hat    ny   estimate of the constant part of theta (4.17)
%   mu2         ny   output of the low-pass filter C(s), i.e. the adaptive
%                    control actually applied (4.23). It is a STATE, not an
%                    algebraic function, because the filter is what decouples
%                    fast estimation from smooth control.
%
% WHY alpha AND beta RATHER THAN theta DIRECTLY.  theta is a nonlinear function
% of (eta, t) and cannot be identified pointwise. Eq (4.17) says that at each
% instant SOME pair (alpha, beta) reproduces it as alpha*||eta|| + beta, so the
% adaptation estimates that pair. The split matters: beta alone could fit theta
% at any single instant, but the ||eta||-proportional term is what lets the
% estimate stay valid as eta changes, which is exactly the regime the
% controller drives the robot through.
%
% -------------------------------------------------------- reset at footstrike
% eta jumps discontinuously at impact -- that is what the impact map does. The
% predictor does not know this, so eta_tilde would register the jump as an
% enormous uncertainty appearing in zero time, and the adaptation (gain 1e4)
% would react to a modelling artifact rather than to model error.
%
% 'reset' therefore re-seeds eta_hat = eta+ so that eta_tilde starts each step
% at zero, under p.l1.reset_predictor.
%
% alpha_hat and beta_hat are ALWAYS carried across the impact, and that is
% deliberate, not an oversight: they describe a property of the robot -- how
% wrong its model is -- and footstrike does not change the robot. Throwing them
% away every step would restart the estimation from scratch at ~3 Hz and the
% controller would never accumulate anything. mu2 is likewise carried, since
% the filter is a physical part of the controller and has no reason to be
% discontinuous.
%
% See also CH4_CTRL_L1, CH4_ODE_RHS, CH4_STEP.

ny = p.ny;

switch lower(action)

    case 'dim'
        out = 5 * ny;

    case 'init'
        eta0 = varargin{1};
        % eta_hat starts ON the real eta (zero prediction error), estimates
        % start at zero (no uncertainty assumed until measured), filter starts
        % relaxed. This is the only initialization consistent with "the
        % controller has learned nothing yet".
        out = [eta0(:); zeros(ny,1); zeros(ny,1); zeros(ny,1)];

    case 'unpack'
        xi = varargin{1};
        out = struct('eta_hat',   xi(1:2*ny), ...
                     'alpha_hat', xi(2*ny+1 : 3*ny), ...
                     'beta_hat',  xi(3*ny+1 : 4*ny), ...
                     'mu2',       xi(4*ny+1 : 5*ny));

    case 'pack'
        s = varargin{1};
        out = [s.eta_hat(:); s.alpha_hat(:); s.beta_hat(:); s.mu2(:)];

    case 'reset'
        xi        = varargin{1};
        eta_plus  = varargin{2};
        out       = xi;
        if p.l1.reset_predictor
            out(1:2*ny) = eta_plus(:);
        end

    otherwise
        error('ch4_l1_state:action', ...
              'Unknown action "%s" (expected dim|init|unpack|pack|reset).', ...
              action);
end

end
