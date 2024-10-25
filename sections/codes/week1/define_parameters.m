function p = define_parameters()

% This function defines the parameters needed for the capital_accumulation.m script

%% Economic Parameters

    % Depreciation rate
    p.delta = 0.05;

    % Exogenous constant investment
    p.Inv = 5;

%% Boundary Conditions

    % Initial capital
    p.K0 = 10;
    
%% Grid Paramters
    
    % The minimum time value in the time domain (start of the simulation)
    p.tmin = 0;
    
    % The maximum time value in the time domain (end of the simulation)
    p.tmax = 300;

    % The number of time steps
    p.I = 100; 

end