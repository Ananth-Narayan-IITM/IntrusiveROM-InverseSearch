function [Urec,Vrec,err_u,err_v,elapsedROM] = solveROMv4(u,v,uBC,vBC,D2,a,b,alpha,beta,Nt,dt,Fu,Fv,x)
    % Without BC lifting gave error of 1.4 and 2% and from plot right BC
    % was not satisified
    % [Uu,Su,~] = svd(u,'econ');
    % [Uv,Sv,~] = svd(v,'econ');
    % BC lifting
    uFluct = u-uBC(x);
    vFluct = v-vBC(x);
    [Uu,Su,~] = svd(uFluct,'econ');
    [Uv,Sv,~] = svd(vFluct,'econ');
    cumEu = cumsum(diag(Su).^2) / sum(diag(Su).^2);
    cumEv = cumsum(diag(Sv).^2) / sum(diag(Sv).^2);
    ru= find(cumEu>= 0.999, 1);
    rv= find(cumEv>= 0.999, 1);
    Uu = Uu(:,1:ru);
    Uv = Uv(:,1:rv);
    % +(Uu'*(alpha.*Uv))
    % +(Uv'*(beta.*Uu))
    AHat = [Uu'*((a.*(D2*Uu))-alpha.*Uu), Uu'*(alpha.*Uv);
        Uv'*(beta.*Uu), Uv'*((b.*(D2*Uv))-beta.*Uv)];
    
    uHat = zeros(ru,Nt);
    vHat = zeros(rv,Nt);
    uHat(:,1) = Uu'*uFluct(:,1);
    vHat(:,1) = Uv'*vFluct(:,1); 
    uvDel = vBC(x)-uBC(x);
    vuDel = uBC(x)-vBC(x);
    tic;
    for n = 2:Nt
        y = [uHat(:,n-1);vHat(:,n-1)];
    
        k1 = AHat * y+[Uu' * (Fu(x,(n-1)*dt)+alpha*uvDel);Uv'*(Fv(x,(n-1)*dt)+beta*vuDel)];
        k2 = AHat * (y + 0.5*dt*k1)+[Uu'*(Fu(x,(n-1)*dt+0.5*dt)+alpha*uvDel);Uv'*(Fv(x,(n-1)*dt+0.5*dt)+beta*vuDel)];
        k3 = AHat * (y + 0.5*dt*k2)+[Uu'*(Fu(x,(n-1)*dt+0.5*dt)+alpha*uvDel);Uv'*(Fv(x,(n-1)*dt+0.5*dt)+beta*vuDel)];
        k4 = AHat * (y + dt*k3)+[Uu'*(Fu(x,(n-1)*dt+dt)+alpha*uvDel);Uv'*(Fv(x,(n-1)*dt+dt)+beta*vuDel)];
        
        yn = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
        uHat(:,n) = yn(1:ru,1);
        vHat(:,n) = yn(ru+1:end,1);
    end
    elapsedROM = toc;
    % Urec = Uu*uHat; Vrec = Uv*vHat;
    % Adapting BC
    Urec = uBC(x)+ Uu*uHat;
    Vrec = vBC(x)+ Uv*vHat;
    err_u = norm(Urec(:,end) - u(:,end), 'fro')*100/norm(u,'fro');
    err_v = norm(Vrec(:,end) - v(:,end), 'fro')*100/norm(v,'fro');
end