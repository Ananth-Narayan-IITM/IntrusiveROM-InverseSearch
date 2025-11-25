function [T,D2,elapsedFOM] = solveFOMv3(TL,TR,Nx,a,t,dt)
    Tint = (TL+TR)/2; x = linspace(0,1,Nx); dx = x(2)-x(1);
    Nt = length(0:dt:t);
    D2 = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
    D2(1,:) = 0; D2(end,:) = 0;
    D2(1,1) = 1; D2(end,end) = 1;
    T = Tint+zeros(Nx,Nt);
    % Full order Model (Nx x Nt)
    tic;
    for n = 2:Nt
        k1 = a * (D2 * T(:,n-1));
        k2 = a * (D2 * (T(:,n-1) + 0.5*dt*k1));
        k3 = a * (D2 * (T(:,n-1) + 0.5*dt*k2));
        k4 = a * (D2 * (T(:,n-1) + dt*k3));
        T(:,n) = T(:,n-1) + dt*(k1 + 2*k2 + 2*k3 + k4) / 6;
    
        % enforce BC after update
        T(1,n) = TL;
        T(Nx,n) = TR;
    end
    elapsedFOM = toc;
end
