function [J,elapsedFOM] = goldenLossFOMv3(a, TObs, TL, TR, Nx, t, dt)
    [T,~,elapsedFOM] = solveFOMv3(TL, TR, Nx, a, t, dt);
    J = 0.5 * norm(T(:,end) - TObs(:,end))^2;
end