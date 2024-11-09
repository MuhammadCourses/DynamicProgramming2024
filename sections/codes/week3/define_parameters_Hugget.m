function p = define_parameters()

% This function defines the parameters needed for the HJB_ramsey_implicit.m script

%% Economic Parameters

    % Discount rate
    p.rho = 0.05;

    % Relative risk aversion coefficient
    p.gamma = 2;

    % Income (Employed and Unemployed)
    p.ze = 0.2;
    p.zu = 0.1; 
    p.z  = [p.ze, p.zu];
    % Hugget transition rates 
    p.lamba_e=0.03;
    p.lamba_u=0.02;
    p.la=[p.lamba_e,p.lamba_u];

    % interest rate 
    p.r = 0.03;

%%    Economic Functions
    
    % Utility function
    p.u = @(c) (c.^(1-p.gamma))./(1-p.gamma);

    % Marginal utility function
    p.mu = @(c) c.^(-p.gamma);

    % Inverse of marginal utility
    % Note: FOC: mu(c) = dv(a) -> c = inv_mu(dv)
    p.inv_mu = @(dv) dv.^(-1/p.gamma);

    % Production function
    p.f = @(k) p.A * k.^p.alpha;

%% Grid Paramters

    % The number of grid points
    p.I = 500;

    % Borrowing Constraint
    p.a_min = -0.02;
    p.a_max = 2;

%% Tuning Parameters
    
    % The maximum number of iteration
    p.maxit = 1000;

    % Convergence criterion, and
    p.tol = 1e-8;

    % Step size (Delta)
    p.Delta = 1000;

end