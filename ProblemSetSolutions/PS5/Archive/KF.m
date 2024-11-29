
function g = KF(s, param)
    % Input: s is the savings policy function (I x 2 matrix)
    % Output: g is stationary distribution (I x 2 matrix)
    
    % Grid parameters
    I = param.I;
    dk = (param.kmax - param.kmin)/(I-1);
    
    % 1. Build transition matrix from savings policy
    A = zeros(2*I, 2*I);
    
    % Loop over states
    for iz = 1:2
        % Index ranges for each state
        ix = (iz-1)*I + 1 : iz*I;
        
        % Forward and backward differences for upwind scheme
        sf = max(s(:,iz), 0);  % Forward drift
        sb = min(s(:,iz), 0);  % Backward drift
        
        % Add drift terms
        for i = 2:I-1
            A(ix(i), ix(i-1)) = -sb(i)/dk;
            A(ix(i), ix(i))   = sf(i)/dk - sb(i)/dk;
            A(ix(i), ix(i+1)) = -sf(i)/dk;
        end
    end
    
    % 2. Add Poisson transition terms
    Az = [-speye(I)*param.la1,  speye(I)*param.la1;
           speye(I)*param.la2, -speye(I)*param.la2];
    A = A + Az;
    
    % 3. Solve for stationary distribution
    AT = A';
    AT(1,:) = 1;  % Normalize probability
    b = zeros(2*I,1);
    b(1) = 1;
    
    % Solve system
    gg = AT \ b;
    
    % Reshape and normalize
    g = reshape(gg, I, 2);
    g = g / (sum(sum(g)) * dk);
    
    % Check if distribution integrates to 1
    if abs(sum(sum(g))*dk - 1) > 1e-5
        warning('Distribution not normalized to 1')
    end
end