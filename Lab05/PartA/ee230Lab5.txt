clear all; close all; clc; 

f = 1:100:1e5;
w = 2*pi*f;
R1 = 330; 
R2 = 1000; 
R3 = 2200; 
L = 47e-3;
Rp = (R1*R2)/(R1+R2);

H = -(R3/R2) ./ ( 1 + 1j*w*L/Rp );
mag_H = abs(H);
phase_H = angle(H)*180/pi; 

semilogx(f, mag_H), grid on
xlabel('Frequency (Hz)'), ylabel('|H(j\omega)|')
title('Magnitude of H(f)')

figure;
semilogx(f, phase_H, 'LineWidth', 1.5);
grid on; hold on;
xline(fp, '--r', sprintf('f_p = %.1f Hz', fp));
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('Phase Response of H(f)');
legend('Phase', 'f_p', 'Location', 'best');