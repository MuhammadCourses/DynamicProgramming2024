function p = define_parameters()

    % This function defines the parameters needed for the Aiyagari model.
    
    %% Economic Parameters
        
        % Relative risk aversion
        p.gamma = 2;
    
        % Discount rate
        p.rho = 0.05;

        % Depreciation rate
        p.delta = 0.05;
        
        % % Initial interest rate
        p.r = 0.025;    
        % capital share
        p.alpha = 1/3;
    
        % exogenous productivity
        p.A = 0.1;
        % Income process
        p.z_u = 1;
        p.z_e = 2;
        p.zz = [p.z_u, p.z_e];
    
        % Probability density   
        p.lambda_u = 1/3;
        p.lambda_e = 1/3;
        p.lambda = [p.lambda_u, p.lambda_e];
        p.L = p.lambda_u/(p.lambda_u + p.lambda_e) * p.z_u + p.lambda_e/(p.lambda_u + p.lambda_e) * p.z_e;
    %% Economic Functions
        
        % Utility funtion
        p.u = @(c) c.^(1-p.gamma)/(1-p.gamma);
    
        % Marginal utility function
        p.mu = @(c) c.^(-p.gamma);
    
        % FOC: mu(c)=dV -> c=inv_mu(dV)
        p.inv_mu = @(dV) dV.^(-1/p.gamma);
    
        % Production function
        p.f = @(x) p.A*x(1)^p.alpha*x(2)^(1-p.alpha);
        p.fk = @(x) p.alpha*p.A*x(1)^(p.alpha-1)*x(2)^(1-p.alpha);
        p.fl = @(x) (1-p.alpha)*p.A*x(1)^p.alpha*x(2)^(-p.alpha);
    
    %% Grid Parmaters
    
        p.kmin = 0;
        p.kmax = 20;
    
        % % The number of grid points
        % p.I = 1000;
    
        % The level of sparse grid
        p.l = 10;
    
    %% Tuning parameters
    
        % Step size: can be arbitrarily large in implicit method
        p.Delta = 1000;
    
        % The maximum number of value function iterations
        p.maxit = 100;
    
        % Tolerance for value function iterations
        p.tol = 10^(-6);
    
        %% TUNING PARAMETERS FOR INTEREST RATE ITERATION
        
        % The maximum number of interest rate iterations
        p.Nr = 40;
    
        % Tolerance for interest rate iterations
        p.tol_S = 10^(-6);
    
    end