%% Initialization of MATLAB
clear variables; close all; clc;
AR = 1.618;
width = 15; height = width/AR;
%% Full Order Model
% Constants
a = 0.1; b = 1; alpha = 0.15; beta = 1; U = 0.25; V = 0.15;
% Boundary Conditions
uL = 0; uR = 1; vL = 0.1; vR = 0.2;
% Solver settings (common)
Nx = 150; x = linspace(0,1,Nx); dx = x(2)-x(1); x = x';
t = 0.5; dt = 0.8*dx^2/(2*max(a,b)); Nt = length(0:dt:t);
% t = 0.5; Nt = 50000; D2 = linspace(0,0.5,Nt); dt = D2(2)-D2(1);
D2 = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
D2(1,:) = 0; D2(end,:) = 0;
D2(1,1) = 1; D2(end,end) = 1;
% Initial Conditions
uI = @(x) uR*x + U*sin(pi*x) - uL*(x - 1);
vI = @(x) vR*x + V*sin(pi*x) - vL*(x - 1);
% Method of manufacturing solutions Refer, MMS_v2
Fu = @(x,t) alpha*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x)) - 2*U*pi*sin(2*pi*t)*sin(pi*x) + U*a*pi^2*cos(2*pi*t)*sin(pi*x);
Fv = @ (x,t) V*b*pi^2*cos(2*pi*t)*sin(pi*x) - 2*V*pi*sin(2*pi*t)*sin(pi*x) - beta*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x));
uVal = uR*x - uL*(x - 1) + U*cos(2*pi*t)*sin(pi*x);
vVal = vR*x - vL*(x - 1) + V*cos(2*pi*t)*sin(pi*x);
% Solving FOM for Validation (Perfect)
[uFOM,vFOM,elapsedFOM] = solveFOMv4(Nx,Nt,dt,uI,vI,uL, uR, vL, vR, D2, alpha, beta, a, b, Fu, Fv);
errFOM = (norm([uFOM(:,end),vFOM(:,end)],'fro')-norm([uVal,vVal],'fro'))*100/norm([uVal,vVal],'fro');
%% Reduced order model
uBC = @(x)(uL*(1 - x) + uR.*x);
vBC = @(x)(vL*(1 - x) + vR.*x);
[uROM,vROM,~,~,elapsedROM] = solveROMv4(uFOM,vFOM,uBC,vBC,D2,a,b,alpha,beta,Nt,dt,Fu,Fv,x);
errROM = (norm([uROM,vROM],'fro')-norm([uFOM,vFOM],'fro'))*100/norm([uFOM,vFOM],'fro');
%% Plotting
% R4v1
% Validation and corresponding plot of FOM, ROM
% Solver common settings
Nx = 150; Lx = 1; x = linspace(0, Lx, Nx)'; dx = x(2)-x(1);
t = 0.5;
D2 = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
D2(1,:) = 0; D2(end,:) = 0; D2(1,1) = 1; D2(end,end) = 1;

% Boundary conditions
uL = 0; uR = 1; vL = 0.1; vR = 0.2;

% 3 sets of (a,b,U,V)
paramSets = [
    0.05  0.1  0.25  0.15;  % Case 1
    0.2   0.5  0.30  0.10;  % Case 2
    1.0   0.8  0.10  0.25]; % Case 3

% Predefine figure
figure('Units','normalized','Position',[0.1 0.1 0.7 0.75]);
tiledlayout(3,2, 'TileSpacing','compact', 'Padding','compact');

for k = 1:3
    a = paramSets(k,1); b = paramSets(k,2);
    U = paramSets(k,3); V = paramSets(k,4);
    alpha = 0.15; beta = 1;

    % Time step and total steps
    dt = 0.8*dx^2/(2*max(a,b)); 
    Nt = length(0:dt:t);

    % Manufactured solutions
    uI = @(x) uR*x + U*sin(pi*x) - uL*(x - 1);
    vI = @(x) vR*x + V*sin(pi*x) - vL*(x - 1);
    Fu = @(x,t) alpha*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x)) ...
                - 2*U*pi*sin(2*pi*t)*sin(pi*x) + U*a*pi^2*cos(2*pi*t)*sin(pi*x);
    Fv = @(x,t) V*b*pi^2*cos(2*pi*t)*sin(pi*x) - 2*V*pi*sin(2*pi*t)*sin(pi*x) ...
                - beta*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x));

    % Analytical
    uVal = uR*x - uL*(x - 1) + U*cos(2*pi*t)*sin(pi*x);
    vVal = vR*x - vL*(x - 1) + V*cos(2*pi*t)*sin(pi*x);

    % --- FOM ---
    [uFOM,vFOM] = solveFOMv4(Nx,Nt,dt,uI,vI,uL,uR,vL,vR,D2,alpha,beta,a,b,Fu,Fv);

    % --- ROM ---
    uBC = @(x)(uL*(1 - x) + uR.*x);
    vBC = @(x)(vL*(1 - x) + vR.*x);
    [uROM,vROM,~,~] = solveROMv4(uFOM,vFOM,uBC,vBC,D2,a,b,alpha,beta,Nt,dt,Fu,Fv,x);

    % --- Compute errors ---
    errFOM(k) = norm([uFOM(:,end),vFOM(:,end)] - [uVal,vVal],'fro')/norm([uVal,vVal],'fro')*100;
    errROM(k) = norm([uROM(:,end),vROM(:,end)] - [uVal, vVal],'fro')/norm([uVal, vVal],'fro')*100;

    % --- Plot u ---
    nexttile;
    ax = gca;
    ax.TitleFontSizeMultiplier = 1;
    ax.FontSize = 12;
    plot(x, uFOM(:,end), 'r-', 'LineWidth', 1.3, 'DisplayName', 'FOM');
    hold on;
    plot(x, uROM(:,end), 'b--', 'LineWidth', 1.3, 'DisplayName', 'ROM');
    plot(x(1:7:end), uVal(1:7:end), 'k*', 'DisplayName', 'Analytical');
    title({sprintf('$a=%.2f$, $b=%.2f$, $U=%.2f$, $V=%.2f$',a,b,U,V)}, 'Interpreter','latex');
    xlabel('$x$', 'Interpreter','latex');
    ylabel('$u(x)$', 'Interpreter','latex');
    % title(sprintf('Case %d: a=%.2f, b=%.2f, U=%.2f, V=%.2f',k,a,b,U,V), 'Interpreter','latex');
    % legend('Location','north');
    grid on;

    % --- Plot v ---
    nexttile;
    ax = gca;
    ax.TitleFontSizeMultiplier = 1;
    ax.FontSize = 12;
    plot(x, vFOM(:,end), 'r-', 'LineWidth', 1.3, 'DisplayName', 'FOM');
    hold on;
    plot(x, vROM(:,end), 'b--', 'LineWidth', 1.3, 'DisplayName', 'ROM');
    plot(x(1:7:end), vVal(1:7:end), 'k*', 'DisplayName', 'Analytical');
    title({sprintf('$a=%.2f$, $b=%.2f$, $U=%.2f$, $V=%.2f$',a,b,U,V)}, 'Interpreter','latex');
    xlabel('$x$', 'Interpreter','latex');
    ylabel('$v(x)$', 'Interpreter','latex');
    % legend('Location','north');
    ax = gca;
    ax.YAxisLocation = 'right';  % Move y-axis to right side
    grid on;

end
lg = legend({'FOM','ROM','Analytical'}, ...
            'Interpreter','latex', ...
            'FontSize',10, ...
            'Orientation','horizontal', ...
            'Box','off');
lg.Layout.Tile = 'south';
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R4f1.pdf', '-dpdf', '-r300');

% R4f2,3
% Performance plot of golden search
k = 10;
aCaptureROM = zeros(1,k);
elapsedFOMInv = zeros(1,k);
aActual = aCaptureROM;
bActual = aCaptureROM;
aCaptureFOM = aCaptureROM;
bCaptureROM = aCaptureROM; bCaptureFOM = aCaptureROM;
elapsedROMInv = aCaptureROM;


for k = 1:k
    a = 0.1+rand*0.9; b = 0.1+rand*1.9;
    % Constants
    alpha = 0.15; beta = 1; U = 0.25; V = 0.15;
    % Boundary Conditions
    uL = 0; uR = 1; vL = 0.1; vR = 0.2;
    % Solver settings (common)
    Nx = 75; x = linspace(0,1,Nx); dx = x(2)-x(1); x = x';
    t = 0.5; dt = 0.8*dx^2/(2*max(a,b)); Nt = length(0:dt:t);
    % t = 0.5; Nt = 50000; D2 = linspace(0,0.5,Nt); dt = D2(2)-D2(1);
    D2 = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
    D2(1,:) = 0; D2(end,:) = 0;
    D2(1,1) = 1; D2(end,end) = 1;
    % Initial Conditions
    uI = @(x) uR*x + U*sin(pi*x) - uL*(x - 1);
    vI = @(x) vR*x + V*sin(pi*x) - vL*(x - 1);
    % Method of manufacturing solutions Refer, MMS_v2
    Fu = @(x,t) alpha*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x)) - 2*U*pi*sin(2*pi*t)*sin(pi*x) + U*a*pi^2*cos(2*pi*t)*sin(pi*x);
    Fv = @ (x,t) V*b*pi^2*cos(2*pi*t)*sin(pi*x) - 2*V*pi*sin(2*pi*t)*sin(pi*x) - beta*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x));
    uVal = uR*x - uL*(x - 1) + U*cos(2*pi*t)*sin(pi*x);
    vVal = vR*x - vL*(x - 1) + V*cos(2*pi*t)*sin(pi*x);
    [uFOM,vFOM] = solveFOMv4(Nx,Nt,dt,uI,vI,uL, uR, vL, vR, D2, alpha, beta, a, b, Fu, Fv);

    uTrue = uFOM(:,end); vTrue = vFOM(:,end); % Given Data
    r = (sqrt(5)-1)/2;
    aInv = rand; bInv = rand; % Initial guess
    tol = 1e-3;
    aPrev = aInv; 
    bPrev = bInv;
    for i = 1:5
        % --- Golden search on 'a' keeping 'b' fixed ---
        a1 = 0; a2 = 1;
        iterA = 0;
        while true
            iterA = iterA + 1;
            aL = a2 - r*(a2 - a1);
            aR = a1 + r*(a2 - a1);
            [JL,time] = goldenLossFOMv4(Nx,Nt,dt,uI,vI,uL,uR,vL,vR,D2,alpha,beta,aL,bInv,Fu,Fv,uTrue,vTrue);
            elapsedFOMInv(k) = time;
            [JR,time] = goldenLossFOMv4(Nx,Nt,dt,uI,vI,uL,uR,vL,vR,D2,alpha,beta,aR,bInv,Fu,Fv,uTrue,vTrue);
            elapsedFOMInv(k) = elapsedFOMInv(k)+time;
            if JL < JR
                a2 = aR;
            else
                a1 = aL;
            end
    
            apInv = (a1 + a2)/2;
            if abs(apInv - aInv) < tol
                aInv = apInv;
                break;
            end
            aInv = apInv;
        end
    
        % --- Golden search on 'b' keeping 'a' fixed ---
        b1 = 0; b2 = 2;
        iterB = 0;
        while true
            iterB = iterB + 1;
            bL = b2 - r*(b2 - b1);
            bR = b1 + r*(b2 - b1);
    
            [JL,time] = goldenLossFOMv4(Nx,Nt,dt,uI,vI,uL,uR,vL,vR,D2,alpha,beta,aInv,bL,Fu,Fv,uTrue,vTrue);
            elapsedFOMInv(k) = elapsedFOMInv(k)+time;
            [JR,time] = goldenLossFOMv4(Nx,Nt,dt,uI,vI,uL,uR,vL,vR,D2,alpha,beta,aInv,bR,Fu,Fv,uTrue,vTrue);
            elapsedFOMInv(k) = elapsedFOMInv(k)+time;
    
            if JL < JR
                b2 = bR;
            else
                b1 = bL;
            end
    
            bpInv = (b1 + b2)/2;
            if abs(bpInv - bInv) < tol
                bInv = bpInv;
                break;
            end
            bInv = bpInv;
    
        end
        if abs(aInv - aPrev) < tol && abs(bInv - bPrev)< tol
            break;
        end
        aPrev = aInv;
        bPrev = bInv;
    end
    aActual(k) = a;
    bActual(k) = b;
    aCaptureFOM(k) = aInv;
    bCaptureFOM(k) = bInv;

    % ROM inverse search
    uBC = @(x)(uL*(1 - x) + uR.*x);
    vBC = @(x)(vL*(1 - x) + vR.*x);
    nSnap = 2;
    uTrue = uFOM(:,round(linspace(1, Nt, nSnap))); 
    vTrue = vFOM(:,round(linspace(1, Nt, nSnap))); % Given Data
    aPrev = aInv; 
    bPrev = bInv;
    for i = 1:5
        % --- Golden search on 'a' keeping 'b' fixed ---
        a1 = 0; a2 = 1;
        iterA = 0;
        while true
            iterA = iterA + 1;
            aL = a2 - r*(a2 - a1);
            aR = a1 + r*(a2 - a1);
            [JL,time] = goldenLossROMv4(uTrue,vTrue,uBC,vBC,D2,aL,bInv,alpha,beta,Nt,dt,Fu,Fv,x);
            elapsedROMInv(k) = time;
            [JR,time] = goldenLossROMv4(uTrue,vTrue,uBC,vBC,D2,aR,bInv,alpha,beta,Nt,dt,Fu,Fv,x);
            elapsedROMInv(k) = elapsedROMInv(k)+time;
            if JL < JR
                a2 = aR;
            else
                a1 = aL;
            end
    
            apInv = (a1 + a2)/2;
            if abs(apInv - aInv) < tol
                aInv = apInv;
                break;
            end
            aInv = apInv;
        end
    
        % --- Golden search on 'b' keeping 'a' fixed ---
        b1 = 0; b2 = 2;
        iterB = 0;
        while true
            iterB = iterB + 1;
            bL = b2 - r*(b2 - b1);
            bR = b1 + r*(b2 - b1);
    
            [JL,time] = goldenLossROMv4(uTrue,vTrue,uBC,vBC,D2,aInv,bL,alpha,beta,Nt,dt,Fu,Fv,x);
            elapsedROMInv(k) = elapsedROMInv(k)+time;            
            [JR,time] = goldenLossROMv4(uTrue,vTrue,uBC,vBC,D2,aInv,bR,alpha,beta,Nt,dt,Fu,Fv,x);
            elapsedROMInv(k) = elapsedROMInv(k)+time;
    
            if JL < JR
                b2 = bR;
            else
                b1 = bL;
            end
            bpInv = (b1 + b2)/2;
            if abs(bpInv - bInv) < tol
                bInv = bpInv;
                break;
            end
            bInv = bpInv;
    
        end
        if abs(aInv - aPrev) < tol && abs(bInv - bPrev)< tol
            break;
        end
        aPrev = aInv;
        bPrev = bInv;
    end

    aCaptureROM(k) = aInv;
    bCaptureROM(k) = bInv;
end

errFOM = zeros(1,k);
errROM = zeros(1,k);

for k = 1:k
    diffFOM = [aActual(k) - aCaptureFOM(k), bActual(k) - bCaptureFOM(k)];
    diffROM = [aActual(k) - aCaptureROM(k), bActual(k) - bCaptureROM(k)];
    denom   = [aActual(k), bActual(k)];
    
    errFOM(k) = 100 * norm(diffFOM, 2) / norm(denom, 2);
    errROM(k) = 100 * norm(diffROM, 2) / norm(denom, 2);
end
idxFOM = errFOM<50; % To filter out large error
figure; hold on; box on;

plot( elapsedFOMInv(idxFOM), 'ro--', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'FOM');
plot( elapsedROMInv(idxFOM), 'bs-', 'LineWidth', 1.2, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'ROM');
ylabel('Elapsed Time $(s)$', 'Interpreter', 'latex');

% Common x-axis and title
xlabel('Case No.', 'Interpreter', 'latex');
title('Elapsed time for inverse problem', 'Interpreter', 'latex');

% Legend and styling
legend('Location', 'east', 'Interpreter', 'latex', 'Box', 'off');

grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R4f2.pdf', '-dpdf', '-r300');

% --- Plotting the errors ---
figure; hold on; box on;

plot(errFOM(idxFOM), 'ro-', 'LineWidth', 1.3, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'FOM');
plot(errROM(idxFOM), 'bs--', 'LineWidth', 1.3, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'ROM');

xlabel('Case No.', 'Interpreter', 'latex');
ylabel('Error (\%)', 'Interpreter', 'latex');
title('Inverse Estimation Error for $(a,b)$', 'Interpreter', 'latex');

legend('Location', 'northeast', 'Interpreter', 'latex', 'Box', 'off');
grid on;

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R4f3.pdf', '-dpdf', '-r300');

% R4f4
% Convext optimization problem: ROM
aList = linspace(0,2,20);
bList = linspace(0,2,20);
JROM = zeros(length(aList), length(bList));

for i = 1:length(aList)
    for j = 1:length(bList)
        [JROM(i,j),~] = goldenLossROMv4(uFOM(:,end),vFOM(:,end),uBC,vBC,D2,aList(i),bList(j),alpha,beta,Nt,dt,Fu,Fv,x);
    end
end

[A,B] = meshgrid(aList,bList);
figure; hold on; box on;
contourf(A, B, JROM); % filled contours colorbar;
xlabel('$a$', 'Interpreter', 'latex');
ylabel('$b$', 'Interpreter', 'latex');
title('Convexity of inverse problem', 'Interpreter', 'latex');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12);
grid on;

% Mark minimum
[minVal, idx] = min(JROM(:));
[aMin, bMin] = deal(A(idx), B(idx));
colorbar;
plot(aMin, bMin, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);
text(aMin+0.02, bMin, 'Minimum', 'Color', 'r', 'FontSize', 11);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R4f4.pdf', '-dpdf', '-r300');

% R4f5
% Performance of FOM and ROM (elapsed time) with Nx
Nx_list = 25:25:600;
nCases = length(Nx_list);

% Arrays for elapsed time
elapsedFOM = zeros(1, nCases);
elapsedROM = zeros(1, nCases);

% Fixed domain
Lx = 1;

% Boundary conditions
uL = 0; uR = 1; 
vL = 0.1; vR = 0.2;

% Fixed parameters for validation (choose one set or average behavior)
a = 0.2; 
b = 0.5;
U = 0.3;
V = 0.1;
alpha = 0.15; 
beta = 1;

% Final time
t = 0.5;

for idx = 1:nCases
    Nx = Nx_list(idx);
    x = linspace(0, Lx, Nx)'; 
    dx = x(2) - x(1);

    % Laplacian matrix
    D2 = spdiags([ones(Nx,1) -2*ones(Nx,1) ones(Nx,1)], [-1 0 1], Nx, Nx) / dx^2;
    D2(1,:) = 0; D2(end,:) = 0;
    D2(1,1) = 1; D2(end,end) = 1;

    % Time stepping
    dt = 0.8 * dx^2 / (2 * max(a, b));
    Nt = length(0:dt:t);

    % Manufactured solutions
    uI = @(x) uR*x + U*sin(pi*x) - uL*(x - 1);
    vI = @(x) vR*x + V*sin(pi*x) - vL*(x - 1);
    Fu = @(x,t) alpha*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + ...
                 U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x)) ...
               -2*U*pi*sin(2*pi*t)*sin(pi*x) + U*a*pi^2*cos(2*pi*t)*sin(pi*x);
    Fv = @(x,t) V*b*pi^2*cos(2*pi*t)*sin(pi*x) - 2*V*pi*sin(2*pi*t)*sin(pi*x) ...
               -beta*(uR*x - vR*x - uL*(x - 1) + vL*(x - 1) + ...
               U*cos(2*pi*t)*sin(pi*x) - V*cos(2*pi*t)*sin(pi*x));

    % --- Run FOM ---
    [uFOM, vFOM, elapsedFOM(idx)] = solveFOMv4(Nx, Nt, dt, uI, vI, ...
                                           uL, uR, vL, vR, D2, ...
                                           alpha, beta, a, b, Fu, Fv);

    % --- Run ROM ---
    uBC = @(x)(uL*(1 - x) + uR.*x);
    vBC = @(x)(vL*(1 - x) + vR.*x);

    [uROM, vROM, ~, ~, elapsedROM(idx)] = solveROMv4(uFOM, vFOM, ...
                                                uBC, vBC, D2, ...
                                                a, b, alpha, beta, ...
                                                Nt, dt, Fu, Fv, x);
end

% Fit FOM: log(T) vs log(Nx)
pFOM = polyfit(log(Nx_list), log(elapsedFOM), 1);
pFOM_slope = pFOM(1);
cFOM = exp(pFOM(2));

% Fit ROM
pROM = polyfit(log(Nx_list), log(elapsedROM), 1);
pROM_slope = pROM(1);
cROM = exp(pROM(2));

fprintf('FOM Complexity: T = %.3e * Nx^{%.2f}\n', cFOM, pFOM_slope);
fprintf('ROM Complexity: T = %.3e * Nx^{%.2f}\n', cROM, pROM_slope);

figure; hold on; box on;

% Plot measured data
plot(Nx_list, elapsedFOM, 'ro-', 'LineWidth', 1.4, ...
    'MarkerFaceColor', 'r', 'DisplayName', 'FOM');
plot(Nx_list, elapsedROM, 'bs-', 'LineWidth', 1.4, ...
    'MarkerFaceColor', 'b', 'DisplayName', 'ROM');

% Plot fitted curves
Nx_fit = linspace(min(Nx_list), max(Nx_list), 500);
plot(Nx_fit, cFOM * Nx_fit.^pFOM_slope, 'r--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('FOM fit: $T=%.1e N_x^{%.2f}$', cFOM, pFOM_slope));
plot(Nx_fit, cROM * Nx_fit.^pROM_slope, 'b--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('ROM fit: $T=%.1e N_x^{%.2f}$', cROM, pROM_slope));

xlabel('$N_x$', 'Interpreter', 'latex');
ylabel('Elapsed Time (s)', 'Interpreter', 'latex');
title('Elapsed Time vs Grid Resolution with Complexity Fits', 'Interpreter', 'latex');

legend('Location','northwest', 'Interpreter', 'latex', 'Box','off');
grid on;

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R4f5.pdf', '-dpdf', '-r300');

