function [J,elapsedROM] = goldenLossROMv3(a, D2,TOrg,Nt,dt,Ts,x)
    [T,~,elapsedROM] = solveROMv3(TOrg,D2,a,Nt,dt,Ts,x);
    J = 0.5 * norm(T(:,end) - TOrg(:,end))^2;
end