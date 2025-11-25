function [J,elapsedFOM] = goldenLossFOMv4(Nx,Nt,dt,uI,vI,uL, uR, vL, vR, D2, alpha, beta, a, b, Fu, Fv,uTrue,vTrue)
    [uFOM,vFOM,elapsedFOM] = solveFOMv4(Nx,Nt,dt,uI,vI,uL, uR, vL, vR, D2, alpha, beta, a, b, Fu, Fv);
    Ju = 0.5 * norm(uTrue(:,end) - uFOM(:,end))^2;
    Jv = 0.5 * norm(vTrue(:,end) - vFOM(:,end))^2;
    J = Ju+Jv;
end