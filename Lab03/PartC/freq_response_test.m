%% Bode fit of a 1st-order stage with possible inversion and gain offset
% Expected CSV columns: Freq, Gain_dB, Phase (Phase in degrees)

clear; clc; close all;

%% === Load data ===
fname = 'EE230PartCTry2.csv';
T = readtable(fname);

% Column names (edit here if your CSV uses different headers)
ff = T.Freq;         % frequency axis (Hz or kHz; auto-checked below)
gain_dB_meas = T.Gain_dB;   % measured magnitude in dB (20*log10(Vout/Vin))
phase_meas_deg = T.Phase;   % measured phase in degrees

%% === Frequency unit sanity check (kHz → Hz if needed) ===
% Heuristic: if max frequency looks like "100" for a sweep that should go to ~100 kHz,
% convert from kHz to Hz.
if max(ff) < 2e3
    ff = ff * 1e3;  % kHz → Hz
end

% Keep only finite points
valid = isfinite(ff) & isfinite(gain_dB_meas) & isfinite(phase_meas_deg);
ff = ff(valid); gain_dB_meas = gain_dB_meas(valid); phase_meas_deg = phase_meas_deg(valid);

% Sort by frequency just in case
[ff, idx] = sort(ff(:));
gain_dB_meas = gain_dB_meas(idx);
phase_meas_deg = phase_meas_deg(idx);

%% === Unwrap measured phase (degrees) for a clean comparison ===
phase_u_deg = unwrap(deg2rad(phase_meas_deg)) * 180/pi;

%% === Quick initial guesses for K and wc ===
% Estimate low-frequency plateau in magnitude (median of lowest decade)
nLF = max(5, round(0.1 * numel(ff)));           % take ~lowest 10% as "LF" region
[~, iLF] = mink(ff, nLF);
LF_mag_dB = median(gain_dB_meas(iLF), 'omitnan');

% Guess |K| from low-frequency mag (dB → linear)
K_abs_guess = 10^(LF_mag_dB/20);
if ~isfinite(K_abs_guess) || K_abs_guess == 0
    K_abs_guess = 1;
end

% Guess sign of K from low-frequency phase (near 0→-90 means +, near -180 means -)
LF_phase_deg = median(phase_u_deg(iLF), 'omitnan');
if ~isfinite(LF_phase_deg), LF_phase_deg = 0; end
if mod(LF_phase_deg + 360, 360) > 90 && mod(LF_phase_deg + 360, 360) < 270
    K_guess = -K_abs_guess;
else
    K_guess = +K_abs_guess;
end

% Guess corner frequency fc by -3 dB point relative to LF plateau (if present)
target_dB = LF_mag_dB - 3;
[~, k] = min(abs(gain_dB_meas - target_dB));
fc_guess = ff(k);
% Fall back if the knee guess is silly
if ~isfinite(fc_guess) || fc_guess <= 0
    fc_guess = max(1e3, min(1e5, median(ff)));  % broad default
end
wc_guess = 2*pi*fc_guess;

%% === Build fitting objective (uses fminsearch; no toolboxes needed) ===
w_meas = 2*pi*ff;

% Weights to balance magnitude (dB) and phase (deg) residuals
% Normalize by the spread in each dataset to keep terms comparable.
mag_scale  = max(5, iqr(gain_dB_meas));      % ~typical dB range
phase_scale = max(45, iqr(phase_u_deg));     % ~typical deg range

obj = @(p) residuals_1pole(p, w_meas, gain_dB_meas, phase_u_deg, mag_scale, phase_scale);

p0 = [K_guess; wc_guess];

% Constrain wc to positive by optimizing its log internally:
% We'll transform inside residuals_1pole (see helper below).
opts = optimset('Display','off', 'MaxFunEvals', 1e5, 'MaxIter', 1e4);
p_fit = fminsearch(@(x) obj([x(1); x(2)]), p0, opts);

K_fit  = p_fit(1);
wc_fit = abs(p_fit(2));      % ensure positive
fc_fit = wc_fit/(2*pi);

%% === Make smooth model curves across the measured range ===
fmin = max(min(ff),  max(1e-3, min(ff))/1.2);
fmax = min(max(ff),  max(ff)*1.2);
fplot = logspace(log10(fmin), log10(fmax), 2000);
wplot = 2*pi*fplot;

% 1st-order LP magnitude and phase with gain K
magH_lin  = abs(K_fit) ./ sqrt(1 + (wplot./wc_fit).^2);
magH_dB   = 20*log10(magH_lin);
phaseLP   = -atan(wplot./wc_fit) * 180/pi;   % base phase
phase_off = (K_fit < 0) * (-180);            % add -180° if inverting
phase_deg = phaseLP + phase_off;

%% === Print results ===
fprintf('Fitted parameters:\n');
fprintf('  K       = %+g (%.2f dB)\n', K_fit, 20*log10(abs(K_fit)));
fprintf('  wc      = %.3e rad/s\n', wc_fit);
fprintf('  fc      = %.3f Hz (%.3f kHz)\n', fc_fit, fc_fit/1e3);

%% === Plots ===
figure(1); clf;
semilogx(fplot, magH_dB, 'LineWidth', 2); hold on;
semilogx(ff, gain_dB_meas, 'x', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on; xlabel('Frequency (Hz)'); ylabel('|H(j\omega)| (dB)');
legend( ...
    sprintf('Model: K=%.2f (%.1f dB), f_c=%.1f kHz', K_fit, 20*log10(abs(K_fit)), fc_fit/1e3), ...
    'Measured', ...
    'Location', 'SouthWest');
set(gca, 'FontSize', 12);
title('Magnitude (dB)');

figure(2); clf;
semilogx(fplot, phase_deg, 'LineWidth', 2); hold on;
semilogx(ff, phase_u_deg, 'x', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on; xlabel('Frequency (Hz)'); ylabel('\angleH(j\omega) (deg)');
legend('Model','Measured (unwrapped)','Location','SouthWest');
set(gca, 'FontSize', 12);
title('Phase (deg)');

%% === Helper (nested) ===
function J = residuals_1pole(p, w, mag_dB_meas, phase_deg_meas, mag_scale, phase_scale)
    % p(1) = K (can be +/-)
    % p(2) = wc (we'll force positive)
    K  = p(1);
    wc = abs(p(2));                 % keep pole positive

    % Model
    mag_dB = 20*log10( abs(K) ./ sqrt(1 + (w./wc).^2) );
    ph     = -atan(w./wc) * 180/pi + (K < 0)*(-180);

    % Residuals with simple balancing between mag and phase
    r_mag = (mag_dB - mag_dB_meas) / mag_scale;
    r_ph  = (ph - phase_deg_meas) / phase_scale;

    % Sum of squares
    J = sum(r_mag.^2 + r_ph.^2);
end
