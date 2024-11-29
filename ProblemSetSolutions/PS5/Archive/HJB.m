function [c, s] = HJB(V, kk, r, w, param)
    % Initialize
    I = param.I;
    dk = (param.kmax - param.kmin)/(I-1);
    num0 = 1e-8;  % numerical 0 for upwind scheme
    
    % Forward and backward differences
    dVf = zeros(I, 2);
    dVb = zeros(I, 2);
    
    for iz = 1:2
        % Interior points
        dVf(1:I-1, iz) = (V(2:I, iz) - V(1:I-1, iz))/dk;
        dVb(2:I, iz) = (V(2:I, iz) - V(1:I-1, iz))/dk;
        
        % Boundary points
        dVf(I, iz) = (w*param.zz(iz) + r*kk(I,iz))/param.gamma;  % state constraint
        dVb(1, iz) = (w*param.zz(iz) + r*kk(1,iz))/param.gamma;  % state constraint
    end
    
    % Consumption from FOC: u'(c) = V'(k)
    cf = param.u1inv(dVf);  % forward difference
    cb = param.u1inv(dVb);  % backward difference
    
    % Savings
    sf = r*kk + w*param.zz - cf;  % forward difference
    sb = r*kk + w*param.zz - cb;  % backward difference
    
    % Upwind scheme
    If = (sf > num0);    % forward difference
    Ib = (sb < -num0);   % backward difference
    I0 = ~If & ~Ib;      % zero drift
    
    % Policy functions
    c = cf.*If + cb.*Ib;
    s = sf.*If + sb.*Ib;
    
    % Make sure consumption is positive
    c = max(c, num0);
    
    % % Verify natural borrowing limit
    % s(kk <= param.kmin,:) = max(0, s(kk <= param.kmin,:));
end