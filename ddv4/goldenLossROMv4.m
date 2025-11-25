function [J,elapsedROM] = goldenLossROMv4(uTrue,vTrue,uBC,vBC,D2,a,b,alpha,beta,Nt,dt,Fu,Fv,x)
    [uROM,vROM,~,~,elapsedROM] = solveROMv4(uTrue,vTrue,uBC,vBC,D2,a,b,alpha,beta,Nt,dt,Fu,Fv,x);
    Ju = 0.5 * norm(uTrue(:,end) - uROM(:,end))^2;
    Jv = 0.5 * norm(vTrue(:,end) - vROM(:,end))^2;
    J = Ju+Jv;
end