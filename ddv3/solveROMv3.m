function [TRecon,err,elapsedROM] = solveROMv3(T,D2,a,Nt,dt,Ts,x) % T = TOrg(:,[1,end]) Initial and final temperatures
    TFluct = T-Ts(x); % Boundary condition lifted
    [U,S,~] = svd(TFluct,"econ");
    cumET = cumsum(diag(S).^2) / sum(diag(S).^2);
    rT= find(cumET>= 0.999, 1);
    % Projection
    U = U(:,1:rT); % Modes
    THat = zeros(rT, Nt);
    THat(:,1) = U'*TFluct(:,1); % Initial Condition
    AHat = a*U'*D2*U;
    tic;
    for n = 2:Nt % Runge Kutta 4th order
        k1 = AHat * THat(:,n-1);
        k2 = AHat * (THat(:,n-1) + 0.5*dt*k1);
        k3 = AHat * (THat(:,n-1) + 0.5*dt*k2);
        k4 = AHat * (THat(:,n-1) + dt*k3);
        THat(:,n) = THat(:,n-1) + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    end
    TRecon = U*THat + Ts(x); % Reconstructed
    err = norm(TRecon(:,end) - T(:,end), 'fro')*100 / norm(T(:,end), 'fro');
    elapsedROM = toc;
end
