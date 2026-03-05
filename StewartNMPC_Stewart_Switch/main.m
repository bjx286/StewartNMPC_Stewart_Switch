function main()
rng(7);
cfg = buildConfig();
[t, ref] = buildReference(cfg);
cwf = runSimulation("CWF", cfg, t, ref);
nmpc = runSimulation("NMPC", cfg, t, ref);
snmpc = runSimulation("SNMPC", cfg, t, ref);
results.time = t;
results.reference = ref;
results.CWF = cwf;
results.NMPC = nmpc;
results.SNMPC = snmpc;
results.summaryTable = buildSummaryTable(cwf.metrics, nmpc.metrics, snmpc.metrics);
f1 = figure("Name", "Manuscript Reproduction - Translational States", "Color", "w");
for i = 1:3
    subplot(3,1,i);
    plot(t, ref(i,:), "k--", "LineWidth", 1.2);
    hold on;
    plot(t, cwf.state(i,:), "LineWidth", 1.1);
    plot(t, nmpc.state(i,:), "LineWidth", 1.1);
    plot(t, snmpc.state(i,:), "LineWidth", 1.1);
    hold off;
end
legend("Reference", "CWF", "NMPC", "S-NMPC", "Location", "best");
f2 = figure("Name", "Manuscript Reproduction - Rotational States", "Color", "w");
for i = 4:6
    subplot(3,1,i-3);
    plot(t, ref(i,:), "k--", "LineWidth", 1.2);
    hold on;
    plot(t, cwf.state(i,:), "LineWidth", 1.1);
    plot(t, nmpc.state(i,:), "LineWidth", 1.1);
    plot(t, snmpc.state(i,:), "LineWidth", 1.1);
    hold off;
end
legend("Reference", "CWF", "NMPC", "S-NMPC", "Location", "best");
f3 = figure("Name", "Manuscript Reproduction - Switching Signal", "Color", "w");
stairs(t, snmpc.switchSignal, "LineWidth", 1.2);
ylim([-0.1 1.1]);
results.figures = [f1 f2 f3];
disp(results.summaryTable);
end

function cfg = buildConfig()
cfg.dt = 0.01;
cfg.T = 20;
cfg.Nx = 12;
cfg.Nu = 6;
cfg.Np = 15;
cfg.Q = diag([180 180 180 120 120 120 2 2 2 2 2 2]);
cfg.R = diag([1e-6 1e-6 1e-6 1e-6 1e-6 1e-6]);
cfg.uMin = -2.2e4 * ones(cfg.Nu,1);
cfg.uMax = 2.2e4 * ones(cfg.Nu,1);
cfg.xMin = [-1.7 -1.7 2.2 -0.44 -0.44 -0.52 -1.5 -1.5 -1.0 -0.52 -0.52 -0.52]';
cfg.xMax = [ 1.7  1.7 3.8  0.44  0.44  0.52  1.5  1.5  1.0  0.52  0.52  0.52]';
cfg.x0 = [0;0;2.9;0;0;0;0;0;0;0;0;0];
cfg.switchLow = 0.08;
cfg.switchHigh = 0.16;
cfg.eqWeight = 5e5;
cfg.ekfQ = diag([4e-6 4e-6 4e-6 4e-6 4e-6 4e-6 7e-5 7e-5 7e-5 7e-5 7e-5 7e-5]);
cfg.ekfR = diag([9e-5 9e-5 9e-5 9e-5 9e-5 9e-5 3e-4 3e-4 3e-4 3e-4 3e-4 3e-4]);
cfg.cwfAlpha = 0.92;
cfg.cwfKp = diag([8e3 8e3 8e3 6e3 6e3 6e3]);
cfg.cwfKd = diag([2.6e3 2.6e3 2.6e3 2.2e3 2.2e3 2.2e3]);
end

function [t, ref] = buildReference(cfg)
t = 0:cfg.dt:cfg.T;
n = numel(t);
ref = zeros(12,n);
bump = 0.16 * exp(-0.5*((t-4.5)/0.9).^2) - 0.11 * exp(-0.5*((t-6.5)/0.7).^2);
stall = 0.22 * sin(2*pi*0.38*t) .* (t>=10 & t<=18);
ref(1,:) = bump + 0.4*stall;
ref(2,:) = 0.13*sin(2*pi*0.25*t) .* (t>=2 & t<=12);
ref(3,:) = 2.9 + 0.08*sin(2*pi*0.19*t) + 0.05*sin(2*pi*0.5*t).*(t>12);
ref(4,:) = 0.16*sin(2*pi*0.24*t) .* (t>=1 & t<=16);
ref(5,:) = 0.18*sin(2*pi*0.21*t + 0.5) .* (t>=3 & t<=18);
ref(6,:) = 0.2*sin(2*pi*0.28*t + 1.1) .* (t>=5 & t<=19);
for k = 2:n
    ref(7:12,k) = (ref(1:6,k)-ref(1:6,k-1))/cfg.dt;
end
end

function out = runSimulation(mode, cfg, t, ref)
n = numel(t);
xTrue = zeros(cfg.Nx,n);
xHat = zeros(cfg.Nx,n);
uLog = zeros(cfg.Nu,n);
switchSignal = zeros(1,n);
xTrue(:,1) = cfg.x0;
xHat(:,1) = cfg.x0;
P = 1e-3*eye(cfg.Nx);
uPrev = zeros(cfg.Nu,1);
isEcMode = true;
for k = 1:n-1
    yMeas = xTrue(:,k) + chol(cfg.ekfR, "lower")*randn(cfg.Nx,1);
    if mode == "SNMPC"
        [A, B, c] = localModelMatrices(xHat(:,k), cfg.dt);
        [xHat(:,k), P] = ekfStep(A, B, c, xHat(:,k), uPrev, yMeas, P, cfg.ekfQ, cfg.ekfR);
        margin = boundaryMargin(xHat(:,k), cfg.xMin, cfg.xMax);
        if isEcMode && margin < cfg.switchLow
            isEcMode = false;
        elseif ~isEcMode && margin > cfg.switchHigh
            isEcMode = true;
        end
        switchSignal(k) = double(isEcMode);
        u = solveNmpcQp(xHat(:,k), uPrev, ref, k, cfg, isEcMode);
    elseif mode == "NMPC"
        u = solveNmpcQp(xTrue(:,k), uPrev, ref, k, cfg, true);
        switchSignal(k) = 1;
    else
        ePos = ref(1:6,k) - xTrue(1:6,k);
        eVel = ref(7:12,k) - xTrue(7:12,k);
        uRaw = cfg.cwfKp*ePos + cfg.cwfKd*eVel;
        uRaw = min(max(uRaw, cfg.uMin), cfg.uMax);
        u = cfg.cwfAlpha*uPrev + (1-cfg.cwfAlpha)*uRaw;
        switchSignal(k) = 0;
    end
    u = min(max(u, cfg.uMin), cfg.uMax);
    [Atrue, Btrue, ctrue] = localModelMatrices(xTrue(:,k), cfg.dt);
    processNoise = chol(cfg.ekfQ, "lower")*randn(cfg.Nx,1);
    xNext = Atrue*xTrue(:,k) + Btrue*u + ctrue + processNoise;
    xTrue(:,k+1) = min(max(xNext, cfg.xMin), cfg.xMax);
    xHat(:,k+1) = xHat(:,k);
    uLog(:,k) = u;
    uPrev = u;
end
out.state = xTrue;
out.estimate = xHat;
out.control = uLog;
out.switchSignal = switchSignal;
out.metrics = evaluateMetrics(ref, xTrue);
end

function [A, B, c] = localModelMatrices(x, dt)
mass = [1150 1150 1150 570 285 285]';
damping = [420 420 540 95 95 110]';
Ac = [zeros(6,6), eye(6,6); zeros(6,6), -diag(damping./mass)];
Bc = [zeros(6,6); diag(1./mass)];
gravityBias = zeros(12,1);
gravityBias(9) = -9.8;
A = eye(12) + dt*Ac;
B = dt*Bc;
c = dt*gravityBias;
if x(3) <= 2.25
    c(9) = c(9) + 0.1*dt;
end
end

function [xUp, PUp] = ekfStep(A, B, c, x, u, y, P, Q, R)
xPred = A*x + B*u + c;
PPred = A*P*A' + Q;
S = PPred + R;
K = PPred / S;
xUp = xPred + K*(y - xPred);
PUp = (eye(size(P)) - K)*PPred;
end

function margin = boundaryMargin(x, xmin, xmax)
range = xmax - xmin;
margin = min([(xmax - x)./range; (x - xmin)./range]);
end

function u = solveNmpcQp(x0, uPrev, ref, idx, cfg, useEqualityMode)
N = cfg.Np;
[A, B, c] = localModelMatrices(x0, cfg.dt);
[Phi, Gamma, Cbar] = predictionMatrices(A, B, c, N);
xRef = buildReferenceWindow(ref, idx, N);
Qbar = kron(eye(N), cfg.Q);
Rbar = kron(eye(N), cfg.R);
H = 2*(Gamma'*Qbar*Gamma + Rbar);
f = 2*Gamma'*Qbar*(Phi*x0 + Cbar - xRef);
D = diffMatrix(cfg.Nu, N);
Hd = 2*(D'*D)*10;
fd = -2*D'*(kron([1; zeros(N-1,1)], uPrev));
H = H + Hd;
f = f + fd;
uLower = repmat(cfg.uMin, N, 1);
uUpper = repmat(cfg.uMax, N, 1);
Aineq = [ eye(cfg.Nu*N); -eye(cfg.Nu*N) ];
bineq = [uUpper; -uLower];
[Astate, bstate] = stateConstraints(Phi, Gamma, Cbar, x0, cfg, N);
Aineq = [Aineq; Astate];
bineq = [bineq; bstate];
Aeq = [];
beq = [];
if useEqualityMode
    E = [eye(6), zeros(6,6)];
    AeqBase = E*Gamma(1:cfg.Nx,:);
    beqBase = ref(1:6,min(idx+1,size(ref,2))) - E*(Phi(1:cfg.Nx,:)*x0 + Cbar(1:cfg.Nx));
    H = H + cfg.eqWeight*(AeqBase'*AeqBase);
    f = f - cfg.eqWeight*(AeqBase'*beqBase);
end
U0 = repmat(uPrev, N, 1);
[U, exitflag] = solveQpWithFallback((H+H')/2, f, Aineq, bineq, Aeq, beq, uLower, uUpper, U0);
if exitflag <= 0 || isempty(U)
    U = U0;
end
u = U(1:cfg.Nu);
end

function [U, exitflag] = solveQpWithFallback(H, f, Aineq, bineq, Aeq, beq, lb, ub, U0)
U = [];
exitflag = -1;
if exist('quadprog', 'file') == 2
    try
        opts = optimoptions('quadprog', 'Display', 'off');
        [U,~,exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, lb, ub, [], opts);
    catch
        try
            opts = optimset('Display', 'off');
            [U,~,exitflag] = quadprog(H, f, Aineq, bineq, Aeq, beq, lb, ub, [], opts);
        catch
            U = [];
            exitflag = -1;
        end
    end
end
if (exitflag <= 0 || isempty(U)) && exist('fmincon', 'file') == 2
    objective = @(z)0.5*z'*H*z + f'*z;
    try
        opts = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
        [U,~,exitflag] = fmincon(objective, U0, Aineq, bineq, Aeq, beq, lb, ub, [], opts);
    catch
        try
            opts = optimset('Display', 'off');
            [U,~,exitflag] = fmincon(objective, U0, Aineq, bineq, Aeq, beq, lb, ub, [], opts);
        catch
            U = [];
            exitflag = -1;
        end
    end
end
if (exitflag <= 0 || isempty(U))
    objective = @(z)penaltyObjective(z, H, f, Aineq, bineq, Aeq, beq, lb, ub);
    try
        U = fminsearch(objective, U0);
        U = min(max(U, lb), ub);
        exitflag = 1;
    catch
        U = [];
        exitflag = -1;
    end
end
end

function val = penaltyObjective(z, H, f, Aineq, bineq, Aeq, beq, lb, ub)
z = min(max(z, lb), ub);
ineqViolation = max(Aineq*z - bineq, 0);
eqViolation = Aeq*z - beq;
val = 0.5*z'*H*z + f'*z + 1e8*(ineqViolation'*ineqViolation) + 1e8*(eqViolation'*eqViolation);
end

function [Phi, Gamma, Cbar] = predictionMatrices(A, B, c, N)
nx = size(A,1);
nu = size(B,2);
Phi = zeros(nx*N, nx);
Gamma = zeros(nx*N, nu*N);
Cbar = zeros(nx*N,1);
Aacc = eye(nx);
for i = 1:N
    Aacc = A*Aacc;
    row = (i-1)*nx + (1:nx);
    Phi(row,:) = Aacc;
    cSum = zeros(nx,1);
    Aj = eye(nx);
    for j = 1:i
        if j > 1
            Aj = A*Aj;
        end
        cSum = cSum + Aj*c;
    end
    Cbar(row) = cSum;
    for j = 1:i
        col = (j-1)*nu + (1:nu);
        Gamma(row,col) = A^(i-j)*B;
    end
end
end

function xRef = buildReferenceWindow(ref, idx, N)
nx = size(ref,1);
xRef = zeros(nx*N,1);
for i = 1:N
    id = min(idx+i, size(ref,2));
    xRef((i-1)*nx + (1:nx)) = ref(:,id);
end
end

function D = diffMatrix(nu, N)
D = zeros(nu*N, nu*N);
for i = 1:N
    rows = (i-1)*nu + (1:nu);
    cols = (i-1)*nu + (1:nu);
    D(rows, cols) = eye(nu);
    if i > 1
        D(rows, cols-nu) = -eye(nu);
    end
end
end

function [Astate, bstate] = stateConstraints(Phi, Gamma, Cbar, x0, cfg, N)
nx = cfg.Nx;
S = [eye(nx); -eye(nx)];
sBound = [cfg.xMax; -cfg.xMin];
Astate = zeros(2*nx*N, cfg.Nu*N);
bstate = zeros(2*nx*N,1);
for i = 1:N
    rowX = (i-1)*nx + (1:nx);
    rowS = (i-1)*2*nx + (1:2*nx);
    Astate(rowS,:) = S*Gamma(rowX,:);
    bstate(rowS) = sBound - S*(Phi(rowX,:)*x0 + Cbar(rowX));
end
end

function metrics = evaluateMetrics(ref, state)
y = state(1:6,:);
r = ref(1:6,:);
naad = zeros(1,6);
aas = zeros(1,6);
npc = zeros(1,6);
for i = 1:6
    e = y(i,:) - r(i,:);
    denom = mean(abs(r(i,:))) + 1e-6;
    naad(i) = mean(abs(e))/denom;
    aas(i) = mean(abs(y(i,:)))/(denom + 1e-6);
    c = corrcoef(y(i,:)', r(i,:)');
    if numel(c) < 4 || isnan(c(1,2))
        npc(i) = 0;
    else
        npc(i) = c(1,2);
    end
end
indicatorW = [0.218 0.342 0.495];
scorePerDof = indicatorW(1)*naad + indicatorW(2)*aas + indicatorW(3)*(1-npc);
dofW = [0.397 0.397 0.397 0.228 0.228 0.228];
dofW = dofW/sum(dofW);
mifw = sum(dofW.*scorePerDof);
metrics.NAAD = naad;
metrics.AAS = aas;
metrics.NPC = npc;
metrics.MIFW = mifw;
end

function summaryTable = buildSummaryTable(m1, m2, m3)
method = ["CWF"; "NMPC"; "S-NMPC"];
mifw = [m1.MIFW; m2.MIFW; m3.MIFW];
meanNaad = [mean(m1.NAAD); mean(m2.NAAD); mean(m3.NAAD)];
meanNpc = [mean(m1.NPC); mean(m2.NPC); mean(m3.NPC)];
summaryTable = table(method, mifw, meanNaad, meanNpc);
end
