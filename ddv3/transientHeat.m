%% Initialization of MATLAB
clear variables; close all; clc;
AR = 1.618;
width = 16; height = width/AR;
%% Heat conduction- Problem Definition
a = 0.0251; % Thermal Diffusivity
TL = 100; TR = 20; Lx = 1;
Nx = 100; x = linspace(0,Lx,Nx)'; dx = x(2)-x(1);
dt = 0.95*dx^2/(2*a);  
t = 0.15; Nt = length(0:dt:t);
%% Validation
Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
nTerms = 200; 
Bn = zeros(nTerms,1);
for n = 1:nTerms
    Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
end
theta = zeros(Nx,1);
for n = 1:nTerms
    theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
end
TVal = Ts+theta;
%% Full Order Model
[TFOM,D2,elapsedFOM] = solveFOMv3(TL,TR,Nx,a,t,dt);
errFOM = (norm(TFOM(:,end),'fro')-norm(TVal,'fro'))*100/norm(TVal,'fro');
%% Reduced Order Model
Ts = @(x) TL*(Lx - x)/Lx + TR*(x)/Lx;
[TROM,errROM,elapsedROM] = solveROMv3(TFOM,D2,a,Nt,dt,Ts,x);
%% Plotting
% R3f1
% Temperature plot for various a along x
TL = 100; TR = 20; Lx = 1;
Nx = 100; x = linspace(0,Lx,Nx)'; dx = x(2)-x(1);
aVals = [0.05, 0.2, 1.0];
colors = {'r', 'b', 'k'};
figure; hold on;
for i = 1:length(aVals)
    a = aVals(i);
    dt = 0.95*dx^2/(2*a);  t = 0.15; Nt = length(0:dt:t);
    Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
    nTerms = 200; 
    Bn = zeros(nTerms,1);
    for n = 1:nTerms
        Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
    end
    theta = zeros(Nx,1);
    for n = 1:nTerms
        theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
    end
    TVal = Ts+theta;
    [TFOM,D2,~] = solveFOMv3(TL,TR,Nx,a,t,dt);
    errFOM(i) = (norm(TFOM(:,end),'fro')-norm(TVal,'fro'))*100/norm(TVal,'fro');
    Ts = @(x) TL*(Lx - x)/Lx + TR*(x)/Lx;
    [TROM,errROM(i),~] = solveROMv3(TFOM,D2,a,Nt,dt,Ts,x);
    
    % Plot results
    plot(x, TFOM(:,end), 'Color', colors{i}, 'LineWidth', 1.5, ...
    'DisplayName', sprintf('FOM, $\\alpha$ = %.2f', a));
    plot(x, TROM(:,end), '--', 'Color', colors{i}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('ROM, $\\alpha$ = %.2f', a));
    
    plot(x(1:10:end), TVal(1:10:end), 'o', 'Color', colors{i}, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('Analytical, $\\alpha$ = %.2f', a));

end
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex');
legend('NumColumns', 3, ...
       'Location', 'southoutside', ...
       'Interpreter', 'latex', ...
       'Box', 'off');
title(sprintf('Temperature Profiles at t = %.2f s', t), 'Interpreter', 'latex');
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f1.pdf', '-dpdf', '-r300');

% R3f2
% Error plot for FOM, ROM with Analytical data
TL = 100; TR = 20; Lx = 1;
Nx = 100; x = linspace(0,Lx,Nx)'; dx = x(2)-x(1);
aVals = linspace(0.1,1,20);
for i = 1:length(aVals)
    a = aVals(i);
    dt = 0.95*dx^2/(2*a); Nt = length(0:dt:t);
    Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
    nTerms = 200; 
    Bn = zeros(nTerms,1);
    for n = 1:nTerms
        Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
    end
    theta = zeros(Nx,1);
    for n = 1:nTerms
        theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
    end
    TVal = Ts+theta;
    [TFOM,D2] = solveFOMv3(TL,TR,Nx,a,t,dt);
    errFOM(i) = (norm(TFOM(:,end),'fro')-norm(TVal,'fro'))*100/norm(TVal,'fro');
    Ts = @(x) TL*(Lx - x)/Lx + TR*(x)/Lx;
    [TROM,~] = solveROMv3(TFOM,D2,a,Nt,dt,Ts,x);
    errROM(i) = (norm(TROM(:,end),'fro')-norm(TVal,'fro'))*100/norm(TVal,'fro');
end
figure; hold on; box on;

% Use smoother lines connecting points
plot(aVals, errFOM, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'FOM');
plot(aVals, errROM, '--sb', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'ROM');

% Labels and title
xlabel('Thermal diffusivity, $a$ (m$^2$/s)', 'Interpreter', 'latex');
ylabel('Relative error (\%)', 'Interpreter', 'latex');
title(sprintf('Error comparison at $t = %.2f$ s', t), 'Interpreter', 'latex');

% Legend and axes
legend('Location', 'southeast', 'Interpreter', 'latex', 'Box', 'off');
xlim([min(aVals) max(aVals)]);
ylim([min([errFOM, errROM]) 2*max([errFOM, errROM])]);
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f2.pdf', '-dpdf', '-r300');

% R3f3
% Effects of boundary lifting
[TBCNL,err] = solveROMv3(TFOM,D2,a,Nt,dt,@(x) 0*TL*(Lx - x)/Lx + TR*(x)/Lx,x); % Boundary non-lifting
[TBCL,~] = solveROMv3(TFOM,D2,a,Nt,dt,@(x) TL*(Lx - x)/Lx + TR*(x)/Lx,x); % Boundary lifting

Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
nTerms = 200; 
Bn = zeros(nTerms,1);
for n = 1:nTerms
    Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
end
theta = zeros(Nx,1);
for n = 1:nTerms
    theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
end
TVal = Ts+theta;

figure;hold on; box on;
plot(x(1:10:end), TVal(1:10:end), 'k*', 'LineWidth', 1.5, ...
    'DisplayName', 'Analytical');
plot(x, TBCNL(:,end), 'r--', 'LineWidth', 1.5, ...
    'DisplayName', 'No boundary lifting');
plot(x, TBCL(:,end), 'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'With boundary lifting');

% Labels and title
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex');
title(sprintf('Effect of Boundary Lifting for ROM at $t = %.2f$ s', t), 'Interpreter', 'latex');

% Legend and axes
legend('Location', 'default', 'Interpreter', 'latex', 'Box', 'off');
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f3.pdf', '-dpdf', '-r300');

% R3f4
% Performance of FOM and ROM (elapsed time) with Nx
NxVal = round(linspace(50,500,75)); 
for i = 1:length(NxVal)
    a = 0.25;
    Nx = NxVal(i);
    x = linspace(0,Lx,Nx)'; dx = x(2)-x(1);
    dt = 0.95*dx^2/(2*a); Nt = length(0:dt:t);
    Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
    nTerms = 200; 
    Bn = zeros(nTerms,1);
    for n = 1:nTerms
        Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
    end
    theta = zeros(Nx,1);
    for n = 1:nTerms
        theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
    end
    TVal = Ts+theta;
    [TFOM,D2,elapsed] = solveFOMv3(TL,TR,Nx,a,t,dt);
    elapsedFOM(i) = elapsed;
    
    tic;
    [TROM,~,elapsed] = solveROMv3(TFOM,D2,a,Nt,dt,@(x) TL*(Lx - x)/Lx + TR*(x)/Lx,x);
    elapsedROM(i) = elapsed;
end
% --- Fit time-complexity: T = c * Nx^p ---
logNx = log(NxVal);

% Fit FOM
logTFOM = log(elapsedFOM);
coeffFOM = polyfit(logNx, logTFOM, 1);
pFOM = coeffFOM(1);
cFOM = exp(coeffFOM(2));

% Fit ROM
logTROM = log(elapsedROM);
coeffROM = polyfit(logNx, logTROM, 1);
pROM = coeffROM(1);
cROM = exp(coeffROM(2));

% Generate smooth Nx curve for fitted complexity line
NxFine = linspace(min(NxVal), max(NxVal), 300);

TFOM_fit = cFOM * NxFine.^pFOM;
TROM_fit = cROM * NxFine.^pROM;

% --- Plot ---
figure; hold on; box on;

% FOM and ROM raw plots
plot(NxVal, elapsedFOM, '-r', 'LineWidth', 1.5, ...
    'DisplayName', 'FOM');
plot(NxVal, elapsedROM, '-b', 'LineWidth', 1.5, ...
    'DisplayName', 'ROM');

% Fitted complexity curves
plot(NxFine, TFOM_fit, '--r', 'LineWidth', 1.2, ...
    'DisplayName', 'FOM fit');
plot(NxFine, TROM_fit, '--b', 'LineWidth', 1.2, ...
    'DisplayName', 'ROM fit');

% Labels and title
xlabel('Domain discretization count', 'Interpreter', 'latex');
ylabel('Elapsed Time (s)', 'Interpreter', 'latex');
title('Performance comparison between FOM and ROM', 'Interpreter', 'latex');

% Legend
legend('Location', 'northeast', 'Interpreter', 'latex', 'Box', 'off');

% --- Add equations on figure ---
eqFOM = sprintf('FOM: $T = %.3e \\; N_x^{%.2f}$', cFOM, pFOM);
eqROM = sprintf('ROM: $T = %.3e \\; N_x^{%.2f}$', cROM, pROM);

text(0.05, 0.95, {eqFOM, eqROM}, ...
    'Units', 'normalized', ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'VerticalAlignment', 'top');

% Axes formatting
xlim([min(NxVal) max(NxVal)]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);

% Save
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f4.pdf', '-dpdf', '-r300');


% R3f5,6
% Golden search performance and predictibility
TL = 100; TR = 20; Lx = 1;
Nx = 100; x = linspace(0,Lx,Nx)'; dx = x(2)-x(1);
aVals = linspace(0.1,1,10);
TAnalytical = zeros(Nx,length(aVals)); TFOMval = TAnalytical; TROMval = TAnalytical;
aInvFOM = zeros(length(aVals),1);
elapsedInvFOM = aInvFOM;
errInvFOM = aInvFOM;
aInvROM = aInvFOM;
elapsedInvROM = aInvFOM;
errInvROM = aInvROM;
for i = 1:length(aVals)
    a = aVals(i);
    dt = 0.95*dx^2/(2*a);  t = 0.15; Nt = length(0:dt:t);
    Ts = TL*(Lx-x)/Lx+TR*x/Lx; T0 = (TL+TR)/2*ones(Nx,1); 
    nTerms = 200; 
    Bn = zeros(nTerms,1);
    for n = 1:nTerms
        Bn(n,1) = 2/Lx*trapz(x,(T0-Ts).*sin(n*pi*x/Lx));
    end
    theta = zeros(Nx,1);
    for n = 1:nTerms
        theta = theta + Bn(n)*sin(n*pi*x/Lx)*exp(-a*(n*pi/Lx)^2*t);
    end
    TVal = Ts+theta;
    TAnalytical(:,i) = TVal;

    [TFOM,D2] = solveFOMv3(TL,TR,Nx,a,t,dt);
    % Inverse Problem for FOM
    a1 = 0; a2 = 1; % Since a \in [0,1]
    r = (sqrt(5)-1)/2;
    
    % Golden Search technique- FOM
    while abs(a2-a1) > 1e-4
        aL = a2 - r*(a2 - a1);
        aR = a1 + r*(a2 - a1);
        [JL,time] = goldenLossFOMv3(aL, TFOM(:,end), TL, TR, Nx, t, dt);
        elapsedInvFOM(i) = elapsedInvFOM(i)+time;
        [JR,time] = goldenLossFOMv3(aR, TFOM(:,end), TL, TR, Nx, t, dt);
        elapsedInvFOM(i) = elapsedInvFOM(i)+time;
        if JL > JR
            a1 = aL;
        else
            a2 = aR;
        end
        aInvFOM(i) = (a1+a2)/2;
    end
    errInvFOM(i) = abs(aInvFOM(i)-a)*100/a;
    if errInvFOM(i) > 5
        error("Inverse for FOM has been failed");
    end
    Ts = @(x) TL*(Lx - x)/Lx + TR*(x)/Lx;
    tic;
    nSnap = 4;
    a1 = 0; a2 = 1; % Since a \in [0,1]
    while (a2 - a1) > 1e-4
        aL = a2 - r*(a2 - a1);
        aR = a1 + r*(a2 - a1);
        [JL,time] = goldenLossROMv3(aL, D2, TFOM(:,round(linspace(2,Nt,nSnap))), Nt, dt,Ts,x);
        elapsedInvROM(i) = elapsedInvROM(i)+time;
        [JR,time] = goldenLossROMv3(aR, D2, TFOM(:,round(linspace(2,Nt,nSnap))), Nt, dt,Ts,x);
        elapsedInvROM(i) = elapsedInvROM(i)+time;
        if JL > JR
            a1 = aL;
        else
            a2 = aR;
        end
        aInvROM(i) = 0.5*(a1 + a2);
    end
    errInvROM(i) = abs(aInvROM(i)-a)*100/a;
    if errInvROM(i) > 5
        error("Inverse for ROM has been failed");
    end
    [TFOM,D2] = solveFOMv3(TL,TR,Nx,aInvFOM(i),t,dt);
    TFOMval(:,i) = TFOM(:,end);
    [TROM,~] = solveROMv3(TFOM(:,round(linspace(2,Nt,nSnap))),D2,aInvROM(i),Nt,dt,Ts,x);
    TROMval(:,i) = TROM(:,end);
end
figure; hold on; box on;
% Use smoother lines connecting points
plot(aVals, elapsedInvFOM, '-r*', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'FOM');
plot(aVals, elapsedInvROM, '--bx', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'ROM');

% Labels and title
xlabel('Thermal Diffusivity, $\alpha$', 'Interpreter', 'latex');
ylabel('Elapsed Time for Golden search (s)', 'Interpreter', 'latex');
title(sprintf('Performance comparison for inverse problem- FOM, ROM'), 'Interpreter', 'latex');

% Legend and axes
legend('Location', 'southeast', 'Interpreter', 'latex', 'Box', 'off');
xlim([min(aVals) max(aVals)]);
grid on; set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f5.pdf', '-dpdf', '-r300');

figure; hold on; box on;
plot(aVals, errInvFOM, 'ro-', 'LineWidth', 1.2, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'FOM');
plot(aVals, errInvROM, 'bs--', 'LineWidth', 1.2, 'MarkerFaceColor', 'b', ...
    'DisplayName', 'ROM');
ylabel('Inverse Error $(\%)$', 'Interpreter', 'latex');
% Common x-axis and title
xlabel('Actual Thermal Diffusivity, $\alpha$', 'Interpreter', 'latex');
title('Inverse Estimation of Thermal Diffusivity', 'Interpreter', 'latex');
% Legend and styling
legend('Location', 'east', 'Interpreter', 'latex', 'Box', 'off');
xlim([min(aVals) max(aVals)]);
grid on;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1);
set(gcf, 'Units', 'centimeters', 'Position', [2 2 width height]);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperPositionMode', 'manual');
print(gcf, 'Plots/R3f6.pdf', '-dpdf', '-r300');