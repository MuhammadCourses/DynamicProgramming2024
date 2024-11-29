
function [V, policy] = VFI(kk, r, w, param)
    % Initialize
    V = zeros(param.I, param.S);  % Value function
    V_old = V;                    % Previous iteration
    policy = zeros(param.I, 2);   % Policy functions
    diff = inf;                   % Convergence check
    it = 1;                       % Iteration counter
    
    % Safety bounds
    V_min = -1e6;  % Minimum allowed value
    V_max = 1e6;   % Maximum allowed value
    
    % Main iteration loop
    while (diff > param.crit && it <= param.maxit)
        % Store old value function
        V_old = V;
        
        % Solve HJB equation
        [c, s] = HJB(V, kk, r, w, param);
        
        % Input validation
        validateattributes(c, {'numeric'}, {'real', 'finite'});
        validateattributes(s, {'numeric'}, {'real', 'finite'});
        
        % Update value function
        V = max(V_min, min(V_max, V));  % Bound values
        
        % Check for complex values
        if ~isreal(V)
            % Debug information
            fprintf('Complex values detected at iteration %d\n', it);
            fprintf('Max imaginary component: %e\n', max(abs(imag(V(:)))));
            error('Complex values in VFI - check numerical stability');
        end
        
        % Check for NaN/Inf
        if any(isnan(V(:))) || any(isinf(V(:)))
            error('NaN or Inf values detected in value function');
        end
        
        % Update convergence check
        diff = max(abs(V(:) - V_old(:)));
        
        % Update iteration counter
        it = it + 1;
    end
    
    % Check convergence
    if it >= param.maxit
        warning('VFI did not converge within maximum iterations');
    end
    
    % Store policy functions
    policy(:,1) = c(:,1);
    policy(:,2) = s(:,1);
end