function param = define_parameters()

% This function defines the parameters needed for the Huggett_GE.m script

%% Economic Parameters
    
    % Relative risk aversion
    param.gamma = 2;                            % although pset calls it sigma, I will stick with this so that it is easier to use previous code. 

    % Discount rate
    param.rho = 0.05;
    
    % depreciation rate
    param.delta = 0.05;

    % share of capital in production
    param.alpha = 1/3;

    % Total factor productivity
    param.A = 0.1;

    % Earnings parameters:
    param.zz  = [1 ,2];     % idiosyncratic shock(efficency units of labor)
    param.la1 = 1/3;        % poisson rate of shock 1
    param.la2 = 1/3;        % poisson rate of shock 2
    param.L   = param.la1/(param.la1+param.la2) * param.zz(1) + param.la2/(param.la1+param.la2) * param.zz(2);  % average labor supply?

%% Economic Functions
    % Utility function
    param.u     = @(c) c.^(1-param.gamma) / (1-param.gamma); 
    % Marginal utility function
    param.u1    = @(c) c.^(-param.gamma);
    % Inverse of marginal utility function
    param.u1inv = @(c) c.^(-1/param.gamma);

    % Production function
    param.f = @(K) param.A*K.^param.alpha*param.L.^(1-param.alpha);

    % marginal product of capital
    param.fK = @(K) param.alpha*param.A*K.^(param.alpha-1)*param.L.^(1-param.alpha);

%% Grid Parmaters

    param.kmin = 0;
    param.kmax = 20;
    param.min = [param.kmin];
    param.max = [param.kmax];

    % The number of grid points
    param.I = 1000;

    % % WE NO LONGER CONSTRUCT GRID POINTS FOR INTEREST RATES
    % % Grid parameters for interest rate
    % param.rmin = -0.05;
    % param.rmax = 0.04;
    % param.Ir = 20;

%%  Tuning parameters for Value Function and Kolmogorov Forward

    % Step size: can be arbitrarily large in implicit method
    param.Delta = 1000;

    % The maximum number of value function iterations
    param.maxit = 100;

    % Tolerance for value function iterations
    param.crit  = 1e-8;

    param.Delta_KF = 1000;
    param.maxit_KF = 100;
    param.crit_KF  = 1e-8;
    %% TUNING PARAMETERS FOR INTEREST RATE ITERATION
    
    % The maximum number of interest rate iterations
    param.Nr = 40;

    % Tolerance for interest rate iterations
    param.tol_S = 10^(-6);

end