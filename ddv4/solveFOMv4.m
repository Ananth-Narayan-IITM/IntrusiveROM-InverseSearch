function [u,v,elapsedFOM] = solveFOMv4(Nx, Nt,dt, uI,vI,uL, uR, vL, vR, D2, alpha, beta, a, b,Su,Sv)
    u = zeros(Nx,Nt);
    v = zeros(Nx,Nt);
    x = linspace(0,1,Nx);
    u(:,1) = uI(x);
    v(:,1) = vI(x);
    
    u(1,:)=uL; u(end,:)=uR;
    v(1,:)=vL; v(end,:)=vR;
    tic;
    for n = 2:Nt
        Fu1 = a * (D2 * u(:,n-1)) + alpha * (v(:,n-1) - u(:,n-1))+Su(x,(n-1)*dt)';
        Fv1 = b * (D2 * v(:,n-1)) + beta  * (u(:,n-1) - v(:,n-1))+Sv(x,(n-1)*dt)';
    
        u2 = u(:,n-1) + 0.5*dt * Fu1;
        v2 = v(:,n-1) + 0.5*dt * Fv1;
        Fu2 = a * (D2 * u2) + alpha * (v2 - u2)+Su(x,(n-1)*dt)';
        Fv2 = b * (D2 * v2) + beta  * (u2 - v2)+Sv(x,(n-1)*dt)';
    
        u3 = u(:,n-1) + 0.5*dt * Fu2;
        v3 = v(:,n-1) + 0.5*dt * Fv2;
        Fu3 = a * (D2 * u3) + alpha * (v3 - u3)+Su(x,(n-1)*dt)';
        Fv3 = b * (D2 * v3) + beta  * (u3 - v3)+Sv(x,(n-1)*dt)';
    
        u4 = u(:,n-1) + dt * Fu3;
        v4 = v(:,n-1) + dt * Fv3;
        Fu4 = a * (D2 * u4) + alpha * (v4 - u4)+Su(x,(n-1)*dt)';
        Fv4 = b * (D2 * v4) + beta  * (u4 - v4)+Sv(x,(n-1)*dt)';
        
        u(:,n) = u(:,n-1)+dt*(Fu1+2*Fu2+2*Fu3+Fu4)/6;
        v(:,n) = v(:,n-1)+dt*(Fv1+2*Fv2+2*Fv3+Fv4)/6;
    
        u(1,n)=uL; u(end,n)=uR;
        v(1,n)=vL; v(end,n)=vR;
    end
    elapsedFOM = toc;
end