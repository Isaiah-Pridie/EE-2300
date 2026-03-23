clear
clc

% Load the measured data file
table = readtable('EE230PartCTry2.csv');

% Extract the frequency, gain in dB, and phase from the table
ff = table.Freq;
gain_dB = table.Gain_dB;
phase = table.Phase;

% Compute the theoretical frequency response of the circuit.
f = linspace(100,10e4,10e3);    % Frequency vector from 100 Hz to 10 kHz

% Write expression for the magnitude of H 
magH = 10000 ./ sqrt(10000^2 + (2*pi*f).^2);    % Units of V/V
magH_dB = 20*log10(magH);    % Convert to dB

% Write expression for the phase of H 
phH = -atan((2*pi*f)/10000);    % Units of radians
phH_deg = rad2deg(phH);         % Convert to degree



figure(1)
semilogx(f, magH_dB, ff, gain_dB, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Frequency (Hz)');
ylabel('|H(j\omega)| (dB)');
legend('Calculated', 'Measured');
grid;
set(gca, 'fontsize', 12);

figure(2)
semilogx(f, phH_deg, ff, phase, 'x', 'linewidth', 2, 'markersize', 10);
xlabel('Frequency (Hz)');
ylabel('\angleH(j\omega) (deg)');
legend('Calculated', 'Measured');
grid;
set(gca, 'fontsize', 12);
