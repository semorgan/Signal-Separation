%% Setup
clear; close all; clc;

fig_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Background Statistics/figures/';
save_fig = 1;

%% Cross correlation example with Ricker wavelets

num_samples = 10000;
seconds = 5;
Fs = num_samples/seconds;

true_mag1 = 100;
true_mag2 = 40;

t = linspace(0,seconds,num_samples);

title_font_size = 14;
axis_font_size = 16;


N = 1000;
[psi,~] = mexihat(-5,5,N);

shift1 = 50;
shift2 = 1500;
s1 = [zeros(1,shift1) psi zeros(1,num_samples - length(psi)-shift1)];
s2 = [zeros(1,shift2) psi zeros(1,num_samples - length(psi)-shift2)];


[r,lags] = xcorr(s1,s2,'coeff');

%Find largest correlation
idx = find(max(r) == r);
%Negative lag means 1st signal arrived first
time_lag = lags(idx);

actual_time_lag = time_lag/Fs;


figure; set(gcf, 'color', 'w'); 
set(gcf,'Position', [441 543 665 254]);
plot(t,s1,'Linewidth',1.5); hold on;box on; grid on;
plot(t,s2,'Linewidth',1.5);
xlabel('Time (sec)', 'FontSize', 14)
title('Impulsive Signals (Ricker Wavelet)', 'FontSize',title_font_size)
lgd = legend('$\vec{s_1}$','$\vec{s_2}$','Interpreter','Latex');
lgd.FontSize = 16;
if save_fig
    saveas(gcf, [fig_folder_path 'impulsive_signals.png'])
end


figure; set(gcf, 'color', 'w'); 
set(gcf,'Position', [441 543 665 254]);
plot([-flip(t) t(1:end-1)],r,'Linewidth',1.5); box on; grid on;
xlabel('Time Lag (sec)', 'FontSize', 14)
title('Cross Correlation', 'FontSize',title_font_size)
ylim([min(r)-.1, max(r)+.1])

if save_fig
    saveas(gcf, [fig_folder_path 'cross_corr.png'])
end

