clear; clc; close all;

K = 1.0;
tau = 10.0;
L = 1.5;

Kp = 2.0;
Ki = 0.4;

tEnd = 200;
dt = 0.05;
t = (0:dt:tEnd)';

r = ones(size(t));

u_min = 0;
u_max = 100;

padeOrder = 3;

s = tf('s');
G0 = K/(tau*s + 1);
[numD, denD] = pade(L, padeOrder);
D = tf(numD, denD);
G = G0 * D;

C = Kp + Ki/s;
T = feedback(C*G, 1);

fprintf('\n========================================\n');
fprintf('BASELINE SYSTEM PERFORMANCE (NOMINAL)\n');
fprintf('========================================\n');

y0_lin = lsim(T, r, t);

simOpts = struct();
simOpts.dt    = dt;
simOpts.u_min = u_min;
simOpts.u_max = u_max;
simOpts.Tt    = tau/10;
simOpts.Kaw   = 1/simOpts.Tt;

[y0, u0, satTime0] = simulate_closed_loop_discrete(G, Kp, Ki, t, r, simOpts, ...
    'distType', 'none', 'distSig', zeros(size(t)));

info0 = safe_stepinfo(y0, t, 1);
fprintf('Rise Time:          %.3f s\n', info0.RiseTime);
fprintf('Settling Time:      %.3f s\n', info0.SettlingTime);
fprintf('Overshoot:          %.3f %%\n', info0.Overshoot);
fprintf('Peak:               %.3f\n', info0.Peak);
fprintf('Steady-State Error: %.5f\n', abs(1 - y0(end)));
fprintf('Saturation Time:    %.3f s\n', satTime0);

[Gm, Pm, Wcg, Wcp] = margin(C*G);

if isfinite(Gm) && Gm > 0
    Gm_dB = 20*log10(Gm);
else
    Gm_dB = Inf;
end

fprintf('\n--- Stability Margins (Nominal Linear Model) ---\n');
if isfinite(Gm_dB)
    fprintf('Gain Margin:  %.2f dB (at %.3f rad/s)\n', Gm_dB, Wcg);
else
    fprintf('Gain Margin:  Inf dB (at %.3f rad/s)\n', Wcg);
end
fprintf('Phase Margin: %.2f deg (at %.3f rad/s)\n', Pm, Wcp);

if Pm < 30
    fprintf('WARNING: Phase margin < 30 deg (poor robustness)\n');
elseif Pm < 45
    fprintf('CAUTION: Phase margin < 45 deg (marginal robustness)\n');
else
    fprintf('GOOD: Phase margin >= 45 deg (robust nominal margins)\n');
end

figure('Name','Baseline Performance','Position',[80 90 1200 750]);
subplot(2,2,1);
plot(t, y0, 'b-', 'LineWidth', 1.8); hold on;
plot(t, r,  'r--','LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Output');
title('Baseline Closed-Loop Response (Discrete)');
legend('y','r','Location','southeast');

subplot(2,2,2);
margin(C*G); grid on;
title('Bode Plot with Stability Margins');

subplot(2,2,3);
nyquist(C*G); grid on;
title('Nyquist Plot');

subplot(2,2,4);
pzmap(T); grid on;
title('Closed-Loop Pole-Zero Map');

figure('Name','Baseline Control Effort','Position',[120 120 1000 450]);
plot(t, u0, 'LineWidth', 1.6);
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Baseline Control Signal u(t)');
ylim([u_min-10 u_max+10]);

fprintf('\n========================================\n');
fprintf('MONTE CARLO ROBUSTNESS ANALYSIS\n');
fprintf('========================================\n');

N = 50;
rng(42);

dK = 0.30;
dTau = 0.40;
dL = 0.20;

fprintf('Testing %d variants...\n', N);
fprintf('Uncertainty: K±%.0f%%, tau±%.0f%%, L±%.0f%%\n', dK*100, dTau*100, dL*100);

crit = struct();
crit.maxOvershoot = 20;
crit.maxSettling = 60;
crit.maxEss = 0.02;
crit.maxSatTime = 15;
crit.mustBeStable = true;

results = table('Size',[N 14], ...
    'VariableTypes', {'double','double','double','logical','logical', ...
                      'double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'K','tau','L','Stable','Pass', ...
                      'RiseTime','SettlingTime','Overshoot','SteadyStateError', ...
                      'PeakOutput','PeakControl','IAE','SatTime','Pm_deg'});

Yall = NaN(length(t), N);
Uall = NaN(length(t), N);

computeMarginsForVariants = true;

for i = 1:N
    K_i   = K   * (1 + (2*rand - 1)*dK);
    tau_i = tau * (1 + (2*rand - 1)*dTau);
    L_i   = L   * (1 + (2*rand - 1)*dL);

    G0_i = K_i/(tau_i*s + 1);
    [numDi, denDi] = pade(L_i, padeOrder);
    Gi = G0_i * tf(numDi, denDi);

    [y, u, satTime] = simulate_closed_loop_discrete(Gi, Kp, Ki, t, r, simOpts, ...
        'distType', 'none', 'distSig', zeros(size(t)));

    stableFlag = all(isfinite(y)) && all(isfinite(u)) && (max(abs(y)) < 1e6);

    if stableFlag
        info = safe_stepinfo(y, t, 1);
        Ess = abs(1 - y(end));
        PeakY = max(abs(y));
        PeakU = max(abs(u));
        IAE = sum(abs(r - y)) * dt;

        if computeMarginsForVariants
            try
                [~, Pm_i, ~, ~] = margin(C*Gi);
            catch
                Pm_i = NaN;
            end
        else
            Pm_i = NaN;
        end

        passFlag = true;
        if crit.mustBeStable
            passFlag = passFlag && stableFlag;
        end
        passFlag = passFlag && (info.Overshoot <= crit.maxOvershoot);
        passFlag = passFlag && (info.SettlingTime <= crit.maxSettling);
        passFlag = passFlag && (Ess <= crit.maxEss);
        passFlag = passFlag && (satTime <= crit.maxSatTime);
    else
        info = struct('RiseTime',NaN,'SettlingTime',NaN,'Overshoot',NaN,'Peak',NaN);
        Ess = NaN; PeakY = NaN; PeakU = NaN; IAE = NaN; satTime = NaN; Pm_i = NaN;
        passFlag = false;
    end

    results{i,:} = [K_i, tau_i, L_i, stableFlag, passFlag, ...
                    info.RiseTime, info.SettlingTime, info.Overshoot, Ess, ...
                    PeakY, PeakU, IAE, satTime, Pm_i];

    Yall(:,i) = y;
    Uall(:,i) = u;
end

stableRows = results(results.Stable,:);
passRows   = results(results.Pass,:);

stableCount = height(stableRows);
passCount   = height(passRows);

fprintf('\n--- Robustness Summary ---\n');
fprintf('Stable Variants: %d / %d (%.1f%%)\n', stableCount, N, 100*stableCount/N);
fprintf('Pass  Variants:  %d / %d (%.1f%%)\n', passCount, N, 100*passCount/N);

if stableCount > 0
    fprintf('\n--- Statistics (Stable Only) ---\n');
    fprintf('Rise Time:        mean %.2f | std %.2f | max %.2f\n', mean(stableRows.RiseTime,'omitnan'), std(stableRows.RiseTime,'omitnan'), max(stableRows.RiseTime,[],'omitnan'));
    fprintf('Settling Time:    mean %.2f | std %.2f | max %.2f\n', mean(stableRows.SettlingTime,'omitnan'), std(stableRows.SettlingTime,'omitnan'), max(stableRows.SettlingTime,[],'omitnan'));
    fprintf('Overshoot:        mean %.2f | std %.2f | max %.2f\n', mean(stableRows.Overshoot,'omitnan'), std(stableRows.Overshoot,'omitnan'), max(stableRows.Overshoot,[],'omitnan'));
    fprintf('Steady-State Err: mean %.4f | std %.4f | max %.4f\n', mean(stableRows.SteadyStateError,'omitnan'), std(stableRows.SteadyStateError,'omitnan'), max(stableRows.SteadyStateError,[],'omitnan'));
    fprintf('Peak Control:     mean %.2f | std %.2f | max %.2f\n', mean(stableRows.PeakControl,'omitnan'), std(stableRows.PeakControl,'omitnan'), max(stableRows.PeakControl,[],'omitnan'));
    fprintf('Sat Time:         mean %.2f | std %.2f | max %.2f\n', mean(stableRows.SatTime,'omitnan'), std(stableRows.SatTime,'omitnan'), max(stableRows.SatTime,[],'omitnan'));
end

figure('Name','Monte Carlo Robustness','Position',[100 50 1300 850]);

subplot(2,3,1);
if stableCount > 0
    plot(t, Yall(:, results.Stable), 'Color', [0.75 0.75 0.75]); hold on;
end
plot(t, y0, 'b-', 'LineWidth', 2);
plot(t, r,  'r--','LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Output');
title(sprintf('Output Overlay (%d Stable / %d)', stableCount, N));
legend('Stable Variants','Baseline','Setpoint','Location','southeast');
ylim([0 1.5]);

subplot(2,3,2);
if stableCount > 0
    plot(t, Uall(:, results.Stable), 'Color', [0.75 0.75 0.75]); hold on;
end
plot(t, u0, 'b-', 'LineWidth', 2);
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Control Effort Overlay');
ylim([u_min-10 u_max+10]);

subplot(2,3,3);
if stableCount > 0
    histogram(stableRows.RiseTime, 15);
    grid on; xlabel('Rise Time [s]'); ylabel('Count'); title('Rise Time Distribution');
else
    axis off; text(0.1,0.5,'No stable variants','FontSize',12);
end

subplot(2,3,4);
if stableCount > 0
    histogram(stableRows.SettlingTime, 15);
    grid on; xlabel('Settling Time [s]'); ylabel('Count'); title('Settling Time Distribution');
else
    axis off; text(0.1,0.5,'No stable variants','FontSize',12);
end

subplot(2,3,5);
if stableCount > 0
    histogram(stableRows.Overshoot, 15);
    grid on; xlabel('Overshoot [%]'); ylabel('Count'); title('Overshoot Distribution');
else
    axis off; text(0.1,0.5,'No stable variants','FontSize',12);
end

subplot(2,3,6);
if stableCount > 0
    scatter3(stableRows.K, stableRows.tau, stableRows.Overshoot, 30, stableRows.Overshoot, 'filled');
    xlabel('Gain K'); ylabel('Time Constant \tau [s]'); zlabel('Overshoot [%]');
    title('Parameter Sensitivity (Stable Only)');
    grid on; view(45,30); colorbar;
else
    axis off; text(0.1,0.5,'No stable variants','FontSize',12);
end

fprintf('\n========================================\n');
fprintf('DISTURBANCE REJECTION TESTS\n');
fprintf('========================================\n');

band = 0.02;

td = 60;
Ad_meas = 0.25;
d_meas = zeros(size(t));
d_meas(t >= td) = Ad_meas;

[yA, uA, satA] = simulate_closed_loop_discrete(G, Kp, Ki, t, r, simOpts, ...
    'distType', 'measurement', 'distSig', d_meas);

yA_meas = yA + d_meas;
[maxDevA, tRecA] = recovery_metrics(t, yA_meas, 1, td, band);

fprintf('\n--- (A) Measurement Disturbance ---\n');
fprintf('Disturbance: +%.0f%% step on measurement at t=%.1f s\n', Ad_meas*100, td);
fprintf('Max deviation: %.4f (%.2f%%)\n', maxDevA, maxDevA*100);
fprintf('Recovery time: %.3f s (to within %.0f%% band)\n', tRecA, band*100);
fprintf('Saturation time: %.3f s\n', satA);

figure('Name','Disturbance Test A: Measurement','Position',[140 140 1000 600]);
subplot(2,1,1);
plot(t, yA_meas, 'LineWidth',1.6); hold on;
plot(t, r, 'r--','LineWidth',1.2);
xline(td,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('y_{meas}');
title('A) Measurement Disturbance Response');
legend('y_{meas}','r','Location','best');

subplot(2,1,2);
plot(t, uA, 'LineWidth',1.6); hold on;
xline(td,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Control Effort u(t)');
ylim([u_min-10 u_max+10]);

tdB = 60;
Ad_load = 0.25;
d_load = zeros(size(t));
d_load(t >= tdB) = Ad_load;

[yB, uB, satB] = simulate_closed_loop_discrete(G, Kp, Ki, t, r, simOpts, ...
    'distType', 'output', 'distSig', d_load);

[maxDevB, tRecB] = recovery_metrics(t, yB, 1, tdB, band);

fprintf('\n--- (B) True Output/Load Disturbance ---\n');
fprintf('Disturbance: +%.0f%% step added to plant output at t=%.1f s\n', Ad_load*100, tdB);
fprintf('Max deviation: %.4f (%.2f%%)\n', maxDevB, maxDevB*100);
fprintf('Recovery time: %.3f s (to within %.0f%% band)\n', tRecB, band*100);
fprintf('Saturation time: %.3f s\n', satB);

figure('Name','Disturbance Test B: Output/Load','Position',[180 180 1000 600]);
subplot(2,1,1);
plot(t, yB, 'LineWidth',1.6); hold on;
plot(t, r, 'r--','LineWidth',1.2);
xline(tdB,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('y');
title('B) True Output/Load Disturbance Response');
legend('y','r','Location','best');

subplot(2,1,2);
plot(t, uB, 'LineWidth',1.6); hold on;
xline(tdB,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Control Effort u(t)');
ylim([u_min-10 u_max+10]);

tdC = 100;
Ad_in = 0.15;
d_in = zeros(size(t));
d_in(t >= tdC) = Ad_in;

[yC, uC_total, satC] = simulate_closed_loop_discrete(G, Kp, Ki, t, r, simOpts, ...
    'distType', 'input', 'distSig', d_in);

[maxDevC, tRecC] = recovery_metrics(t, yC, 1, tdC, band);

fprintf('\n--- (C) Input/Actuator Disturbance ---\n');
fprintf('Disturbance: +%.0f%% step added at plant input at t=%.1f s\n', Ad_in*100, tdC);
fprintf('Max deviation: %.4f (%.2f%%)\n', maxDevC, maxDevC*100);
fprintf('Recovery time: %.3f s (to within %.0f%% band)\n', tRecC, band*100);
fprintf('Saturation time: %.3f s\n', satC);

figure('Name','Disturbance Test C: Input/Actuator','Position',[220 220 1000 600]);
subplot(2,1,1);
plot(t, yC, 'LineWidth',1.6); hold on;
plot(t, r, 'r--','LineWidth',1.2);
xline(tdC,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('y');
title('C) Input/Actuator Disturbance Response');
legend('y','r','Location','best');

subplot(2,1,2);
plot(t, uC_total, 'LineWidth',1.6); hold on;
plot(t, d_in, ':','LineWidth',1.3);
xline(tdC,'k--','Disturbance');
grid on; xlabel('Time [s]'); ylabel('Signal [%]');
title('Total Input (Saturated)');
legend('u_{total}','d_{in}','Location','best');
ylim([u_min-10 u_max+10]);

fprintf('\n========================================\n');
fprintf('FINAL SUMMARY\n');
fprintf('========================================\n');

summary = table();

summary.Baseline_PhaseMargin_deg = Pm;
summary.Baseline_GainMargin_dB = Gm_dB;
summary.Baseline_RiseTime_s = info0.RiseTime;
summary.Baseline_SettlingTime_s = info0.SettlingTime;
summary.Baseline_Overshoot_pct = info0.Overshoot;
summary.Baseline_Ess = abs(1 - y0(end));
summary.Baseline_SatTime_s = satTime0;

summary.MC_StableRate_pct = 100*stableCount/N;
summary.MC_PassRate_pct = 100*passCount/N;

if stableCount > 0
    summary.MC_WorstOvershoot_pct = max(stableRows.Overshoot,[],'omitnan');
    summary.MC_WorstSettlingTime_s = max(stableRows.SettlingTime,[],'omitnan');
    summary.MC_WorstEss = max(stableRows.SteadyStateError,[],'omitnan');
    summary.MC_WorstSatTime_s = max(stableRows.SatTime,[],'omitnan');
else
    summary.MC_WorstOvershoot_pct = NaN;
    summary.MC_WorstSettlingTime_s = NaN;
    summary.MC_WorstEss = NaN;
    summary.MC_WorstSatTime_s = NaN;
end

summary.MeasDist_MaxDev = maxDevA;
summary.MeasDist_Rec_s = tRecA;

summary.LoadDist_MaxDev = maxDevB;
summary.LoadDist_Rec_s = tRecB;

summary.InDist_MaxDev = maxDevC;
summary.InDist_Rec_s = tRecC;

disp(summary);

doExport = false;
if doExport
    writetable(results, 'monte_carlo_results_v2.csv');
    writetable(summary, 'summary_v2.csv');
end

fprintf('\nANALYSIS COMPLETE\n');

function [y, u_total, satTime] = simulate_closed_loop_discrete(G, Kp, Ki, t, r, simOpts, varargin)
    p = inputParser;
    p.addParameter('distType','none',@(x)ischar(x) || (isstring(x) && isscalar(x)));
    p.addParameter('distSig',zeros(size(t)),@(x)isnumeric(x) && numel(x)==numel(t));
    p.parse(varargin{:});

    distType = lower(char(p.Results.distType));
    d = p.Results.distSig(:);

    dt = simOpts.dt;
    u_min = simOpts.u_min;
    u_max = simOpts.u_max;
    Kaw = simOpts.Kaw;

    sys = ss(G);
    sysd = c2d(sys, dt, 'tustin');

    nx = size(sysd.A,1);
    x = zeros(nx,1);

    y = zeros(size(t));
    u_total = zeros(size(t));

    integE = 0;
    satTime = 0;

    y_prev = 0;

    for k = 1:numel(t)
        if strcmp(distType,'measurement')
            y_meas = y_prev + d(max(k-1,1));
        else
            y_meas = y_prev;
        end

        e = r(k) - y_meas;
        u_raw = Kp*e + Ki*integE;

        if strcmp(distType,'input')
            u_raw_total = u_raw + d(k);
            u_sat = min(max(u_raw_total, u_min), u_max);
            u_total(k) = u_sat;
            integE = integE + dt*( e + Kaw*(u_sat - u_raw_total) );
        else
            u_sat = min(max(u_raw, u_min), u_max);
            u_total(k) = u_sat;
            integE = integE + dt*( e + Kaw*(u_sat - u_raw) );
        end

        if u_total(k) <= u_min+eps || u_total(k) >= u_max-eps
            satTime = satTime + dt;
        end

        x = sysd.A*x + sysd.B*u_total(k);
        y_plant = sysd.C*x + sysd.D*u_total(k);

        if strcmp(distType,'output')
            y(k) = y_plant + d(k);
        else
            y(k) = y_plant;
        end

        y_prev = y(k);

        if ~isfinite(y(k)) || abs(y(k)) > 1e6
            y(k:end) = NaN;
            u_total(k:end) = NaN;
            return;
        end
    end
end

function info = safe_stepinfo(y, t, finalValue)
    if any(~isfinite(y))
        info = struct('RiseTime',NaN,'SettlingTime',NaN,'Overshoot',NaN,'Peak',NaN);
        return;
    end
    try
        si = stepinfo(y, t, finalValue);
        info = struct('RiseTime',si.RiseTime,'SettlingTime',si.SettlingTime, ...
                      'Overshoot',si.Overshoot,'Peak',si.Peak);
    catch
        info = struct('RiseTime',NaN,'SettlingTime',NaN,'Overshoot',NaN,'Peak',NaN);
    end
end

function [maxDev, tRecover] = recovery_metrics(t, y, ref, tDist, band)
    idx = find(t >= tDist, 1, 'first');
    if isempty(idx)
        maxDev = NaN;
        tRecover = NaN;
        return;
    end
    dev = abs(y(idx:end) - ref);
    maxDev = max(dev, [], 'omitnan');

    recIdx = find(dev <= band, 1, 'first');
    if isempty(recIdx)
        tRecover = Inf;
    else
        tRecover = t(idx + recIdx - 1) - tDist;
    end
end