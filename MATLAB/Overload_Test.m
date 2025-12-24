clear; clc; close all;

K = 1.0;
tau = 10.0;
L = 1.5;

Kp = 2.0;
Ki = 0.4;

tEnd = 200;
dt = 0.05;
t = (0:dt:tEnd)';

u_min = 0;
u_max = 100;

padeOrder = 3;

s = tf('s');
G_nom = build_plant_tf(K, tau, L, padeOrder);

C = Kp + Ki/s;
T_nom = feedback(C*G_nom, 1);

r_nom = ones(size(t));

simOpts = struct();
simOpts.dt    = dt;
simOpts.u_min = u_min;
simOpts.u_max = u_max;
simOpts.Tt    = tau/10;
simOpts.Kaw   = 1/simOpts.Tt;

fprintf('\n========================================\n');
fprintf('BASELINE SYSTEM PERFORMANCE (NOMINAL)\n');
fprintf('========================================\n');

[y0, u0, satTime0] = simulate_closed_loop_discrete(G_nom, Kp, Ki, t, r_nom, simOpts, ...
    'distType', 'none', 'distSig', zeros(size(t)));

info0 = safe_stepinfo(y0, t, 1);
fprintf('Rise Time:          %.3f s\n', info0.RiseTime);
fprintf('Settling Time:      %.3f s\n', info0.SettlingTime);
fprintf('Overshoot:          %.3f %%\n', info0.Overshoot);
fprintf('Peak:               %.3f\n', info0.Peak);
fprintf('Steady-State Error: %.5f\n', abs(1 - y0(end)));
fprintf('Saturation Time:    %.3f s\n', satTime0);

[Gm, Pm, Wcg, Wcp] = margin(C*G_nom);
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

figure('Name','Nominal Baseline Summary','Position',[80 90 1200 750]);
subplot(2,2,1);
plot(t, y0, 'b-', 'LineWidth', 1.8); hold on;
plot(t, r_nom, 'r--', 'LineWidth', 1.2);
grid on; xlabel('Time [s]'); ylabel('Output');
title('Baseline Closed-Loop Response (Discrete)');
legend('y','r','Location','southeast');

subplot(2,2,2);
margin(C*G_nom); grid on;
title(sprintf('Bode Plot with Stability Margins (Pm=%.1f deg)', Pm));

subplot(2,2,3);
nyquist(C*G_nom); grid on;
title('Nyquist Plot (Nominal)');

subplot(2,2,4);
pzmap(T_nom); grid on;
title('Closed-Loop Pole-Zero Map (Nominal Linear)');

figure('Name','Nominal Control Signal','Position',[120 120 1000 450]);
plot(t, u0, 'LineWidth', 1.6);
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Baseline Control Signal u(t)');
ylim([u_min-10 u_max+10]);

fprintf('\n========================================\n');
fprintf('OVERLOAD + DEGRADED DYNAMICS TEST\n');
fprintf('========================================\n');

t_ov_start = 60;
t_ov_end   = 140;

r_ov_level = 1.5;
r_overload = r_nom;
r_overload(t >= t_ov_start & t <= t_ov_end) = r_ov_level;

K_ov   = 0.70*K;
tau_ov = 1.60*tau;
L_ov   = 1.20*L;

G_ov = build_plant_tf(K_ov, tau_ov, L_ov, padeOrder);

[y_ov, u_ov, satTime_ov, whichPlant] = simulate_closed_loop_piecewise(G_nom, G_ov, Kp, Ki, t, r_overload, simOpts, ...
    t_ov_start, t_ov_end);

band = 0.02;

idxA = find(t >= t_ov_start, 1, 'first');
idxB = find(t >= t_ov_end,   1, 'first');

if isempty(idxA), idxA = 1; end
if isempty(idxB), idxB = numel(t); end

maxDev_overload = max(abs(y_ov(idxA:idxB) - r_ov_level), [], 'omitnan');

[tRec_after, maxDev_after] = recovery_after_window(t, y_ov, 1, t_ov_end, band);

satTime_overload = sum((u_ov(idxA:idxB) <= u_min+eps) | (u_ov(idxA:idxB) >= u_max-eps)) * dt;

fprintf('Overload window:   [%.1f, %.1f] s\n', t_ov_start, t_ov_end);
fprintf('Setpoint during overload: %.2f\n', r_ov_level);
fprintf('Degraded plant: K=%.2f (from %.2f), tau=%.2f (from %.2f), L=%.2f (from %.2f)\n', ...
    K_ov, K, tau_ov, tau, L_ov, L);

fprintf('\n--- Overload Metrics ---\n');
fprintf('Max deviation DURING overload (vs %.2f): %.4f (%.2f%% of overload SP)\n', ...
    r_ov_level, maxDev_overload, 100*maxDev_overload/r_ov_level);
fprintf('Saturation time DURING overload: %.3f s\n', satTime_overload);
fprintf('Total saturation time (whole test): %.3f s\n', satTime_ov);

fprintf('\n--- Post-Overload Recovery (back to r=1) ---\n');
fprintf('Max deviation AFTER overload: %.4f (%.2f%%)\n', maxDev_after, 100*maxDev_after);
if isfinite(tRec_after)
    fprintf('Recovery time AFTER overload (within %.1f%% band): %.3f s\n', 100*band, tRec_after);
else
    fprintf('Recovery time AFTER overload: Inf (never enters band)\n');
end

figure('Name','Overload Scenario Response','Position',[140 140 1200 700]);
subplot(2,1,1);
plot(t, y_ov, 'b-', 'LineWidth', 1.6); hold on;
plot(t, r_overload, 'r--', 'LineWidth', 1.2);
xline(t_ov_start,'k--','Overload Start');
xline(t_ov_end,'k--','Overload End');
grid on; xlabel('Time [s]'); ylabel('Output');
title('Closed-Loop Output Under Overload + Degraded Dynamics');
legend('y','r','Location','best');

subplot(2,1,2);
plot(t, u_ov, 'LineWidth', 1.6); hold on;
yline(u_min,'k:'); yline(u_max,'k:');
xline(t_ov_start,'k--','Overload Start');
xline(t_ov_end,'k--','Overload End');
grid on; xlabel('Time [s]'); ylabel('u [%]');
title('Control Effort Under Overload (with Saturation Limits)');
ylim([u_min-10 u_max+10]);

figure('Name','Which Plant Was Active','Position',[200 200 1000 250]);
stairs(t, whichPlant, 'LineWidth', 1.6);
grid on; xlabel('Time [s]'); ylabel('Plant ID');
yticks([0 1]); yticklabels({'Nominal','Overload'});
title('Piecewise Plant Switching Indicator');

fprintf('\n========================================\n');
fprintf('MONTE CARLO UNDER OVERLOAD (OPTIONAL)\n');
fprintf('========================================\n');

doMC_overload = true;

if doMC_overload
    N = 50;
    rng(7);

    dK   = 0.30;
    dTau = 0.40;
    dLmc = 0.20;

    resultsOV = table('Size',[N 12], ...
        'VariableTypes', {'double','double','double','logical', ...
                          'double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'K','tau','L','Stable', ...
                          'MaxDevDuring','SatTimeDuring','MaxDevAfter','RecAfter_s', ...
                          'PeakY','PeakU','IAE_total','SatTime_total'});

    for i = 1:N
        K_i   = K   * (1 + (2*rand - 1)*dK);
        tau_i = tau * (1 + (2*rand - 1)*dTau);
        L_i   = L   * (1 + (2*rand - 1)*dLmc);

        G_nom_i = build_plant_tf(K_i, tau_i, L_i, padeOrder);

        K_ov_i   = 0.70*K_i;
        tau_ov_i = 1.60*tau_i;
        L_ov_i   = 1.20*L_i;
        G_ov_i = build_plant_tf(K_ov_i, tau_ov_i, L_ov_i, padeOrder);

        [y_i, u_i, satTot_i] = simulate_closed_loop_piecewise(G_nom_i, G_ov_i, Kp, Ki, t, r_overload, simOpts, ...
            t_ov_start, t_ov_end);

        stableFlag = all(isfinite(y_i)) && all(isfinite(u_i)) && (max(abs(y_i)) < 1e6);

        idxA = find(t >= t_ov_start, 1, 'first');
        idxB = find(t >= t_ov_end,   1, 'first');
        if isempty(idxA), idxA = 1; end
        if isempty(idxB), idxB = numel(t); end

        if stableFlag
            maxDevDuring = max(abs(y_i(idxA:idxB) - r_ov_level), [], 'omitnan');
            satDuring = sum((u_i(idxA:idxB) <= u_min+eps) | (u_i(idxA:idxB) >= u_max-eps)) * dt;

            [recAfter, maxDevAfter_i] = recovery_after_window(t, y_i, 1, t_ov_end, band);

            peakY = max(abs(y_i), [], 'omitnan');
            peakU = max(abs(u_i), [], 'omitnan');

            iaeTotal = sum(abs(r_overload - y_i)) * dt;
        else
            maxDevDuring = NaN; satDuring = NaN; recAfter = NaN; maxDevAfter_i = NaN;
            peakY = NaN; peakU = NaN; iaeTotal = NaN;
        end

        resultsOV{i,:} = [K_i, tau_i, L_i, stableFlag, ...
                          maxDevDuring, satDuring, maxDevAfter_i, recAfter, ...
                          peakY, peakU, iaeTotal, satTot_i];
    end

    stableOV = resultsOV(resultsOV.Stable,:);
    fprintf('MC Overload stable: %d/%d (%.1f%%)\n', height(stableOV), N, 100*height(stableOV)/N);

    if height(stableOV) > 0
        fprintf('MaxDevDuring (mean ± std): %.4f ± %.4f\n', mean(stableOV.MaxDevDuring,'omitnan'), std(stableOV.MaxDevDuring,'omitnan'));
        fprintf('SatTimeDuring (mean ± std): %.3f ± %.3f s\n', mean(stableOV.SatTimeDuring,'omitnan'), std(stableOV.SatTimeDuring,'omitnan'));
        fprintf('RecAfter_s (mean ± std): %.3f ± %.3f s\n', mean(stableOV.RecAfter_s,'omitnan'), std(stableOV.RecAfter_s,'omitnan'));
        fprintf('Worst RecAfter_s: %.3f s\n', max(stableOV.RecAfter_s,[],'omitnan'));
    end

    figure('Name','MC Overload Distributions','Position',[260 120 1200 700]);
    subplot(2,2,1); histogram(stableOV.MaxDevDuring, 15); grid on; xlabel('Max Dev During'); ylabel('Count'); title('Max Deviation During Overload');
    subplot(2,2,2); histogram(stableOV.SatTimeDuring, 15); grid on; xlabel('SatTime During [s]'); ylabel('Count'); title('Saturation Time During Overload');
    subplot(2,2,3); histogram(stableOV.RecAfter_s, 15); grid on; xlabel('Recovery After [s]'); ylabel('Count'); title('Recovery Time After Overload');
    subplot(2,2,4); histogram(stableOV.IAE_total, 15); grid on; xlabel('IAE total'); ylabel('Count'); title('Total IAE (whole test)');
end

fprintf('\n========================================\n');
fprintf('ANALYSIS COMPLETE\n');
fprintf('========================================\n');

function G = build_plant_tf(K, tau, L, padeOrder)
s = tf('s');
G0 = K/(tau*s + 1);
[numD, denD] = pade(L, padeOrder);
D = tf(numD, denD);
G = G0 * D;
end

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

function [y, u_total, satTime, whichPlant] = simulate_closed_loop_piecewise(G_nom, G_ov, Kp, Ki, t, r, simOpts, tStart, tEnd)
dt = simOpts.dt;
u_min = simOpts.u_min;
u_max = simOpts.u_max;
Kaw = simOpts.Kaw;

sysNom = ss(G_nom);
sysdNom = c2d(sysNom, dt, 'tustin');

sysOv = ss(G_ov);
sysdOv = c2d(sysOv, dt, 'tustin');

nxNom = size(sysdNom.A,1);
nxOv  = size(sysdOv.A,1);

xNom = zeros(nxNom,1);
xOv  = zeros(nxOv,1);

y = zeros(size(t));
u_total = zeros(size(t));
whichPlant = zeros(size(t));

integE = 0;
satTime = 0;
y_prev = 0;

for k = 1:numel(t)
    if t(k) >= tStart && t(k) <= tEnd
        active = 1;
    else
        active = 0;
    end
    whichPlant(k) = active;

    e = r(k) - y_prev;
    u_raw = Kp*e + Ki*integE;
    u_sat = min(max(u_raw, u_min), u_max);
    u_total(k) = u_sat;
    integE = integE + dt*( e + Kaw*(u_sat - u_raw) );

    if u_total(k) <= u_min+eps || u_total(k) >= u_max-eps
        satTime = satTime + dt;
    end

    if active == 1
        xOv = sysdOv.A*xOv + sysdOv.B*u_total(k);
        yk = sysdOv.C*xOv + sysdOv.D*u_total(k);
    else
        xNom = sysdNom.A*xNom + sysdNom.B*u_total(k);
        yk = sysdNom.C*xNom + sysdNom.D*u_total(k);
    end

    y(k) = yk;
    y_prev = y(k);

    if ~isfinite(y(k)) || abs(y(k)) > 1e6
        y(k:end) = NaN;
        u_total(k:end) = NaN;
        whichPlant(k:end) = whichPlant(k);
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

function [tRecover, maxDevAfter] = recovery_after_window(t, y, ref, tWindowEnd, band)
idx = find(t >= tWindowEnd, 1, 'first');
if isempty(idx)
    tRecover = NaN;
    maxDevAfter = NaN;
    return;
end
dev = abs(y(idx:end) - ref);
maxDevAfter = max(dev, [], 'omitnan');

recIdx = find(dev <= band, 1, 'first');
if isempty(recIdx)
    tRecover = Inf;
else
    tRecover = t(idx + recIdx - 1) - tWindowEnd;
end
end
