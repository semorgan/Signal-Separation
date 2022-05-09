%% Setup
clear; close all; clc;


ex_name = 'coh_square_hat_shift';
figure_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Array Coherence/figures/';

save_fig = 0;

axis_font_size = 16;
title_font_size = 14;

%% Square and hat functions

a = 0;
b = 1;
A1 = 1;
A2 = 1;



%shifted
Fs = 50000;
num_samples = Fs*(b-a);
t = linspace(a,b,Fs*(b-a));
padding = Fs/2;


s1 = A1*ones(size(t));
s1_nopadding = s1;
s2_temp1 = (2*A2*t(1:length(t)/2)/(b-a))+(2*A2*a)/(a-b);
s2_temp2 =  (2*A2*t(length(t)/2+1:end)/(a-b))+(2*A2*b)/(b-a);
s2_nopadding = [s2_temp1 s2_temp2];
s1 = [zeros(1,padding) s1 zeros(1,padding)];
s2 = [zeros(1,padding) s2_temp1 s2_temp2 zeros(1,padding)];

%Shifting 1st source
speed1 = 343;
%1st shift (15 meters)
shift_val1 = ceil((15/speed1)*Fs);
s1_1 = [zeros(1,shift_val1+padding) s1_nopadding zeros(1,padding - shift_val1)];
%2nd shift (10 meters)
shift_val2 = ceil((10/speed1)*Fs);
s1_2 = [zeros(1,shift_val2+padding) s1_nopadding zeros(1,padding - shift_val2)];
%3rd shift (15 meters)
shift_val3 = ceil((15/speed1)*Fs);
s1_3 = [zeros(1,shift_val3+padding) s1_nopadding zeros(1,padding - shift_val3)];
%4th shift (20 meters)
shift_val4 = ceil((20/speed1)*Fs);
s1_4 = [zeros(1,shift_val4+padding) s1_nopadding zeros(1,padding - shift_val4)];


%Shifting 2nd source
speed2 = 3*10^8;
%1st shift (15 meters)
shift_val1 = ceil((20/speed2)*Fs);
s2_1 = [zeros(1,shift_val1+padding) s2_nopadding zeros(1,padding - shift_val1)];
%2nd shift (10 meters)
shift_val2 = ceil((15/speed2)*Fs);
s2_2 = [zeros(1,shift_val2+padding) s2_nopadding zeros(1,padding - shift_val2)];
%3rd shift (15 meters)
shift_val3 = ceil((10/speed2)*Fs);
s2_3 = [zeros(1,shift_val3+padding) s2_nopadding zeros(1,padding - shift_val3)];
%4th shift (20 meters)
shift_val4 = ceil((15/speed2)*Fs);
s2_4 = [zeros(1,shift_val4+padding) s2_nopadding zeros(1,padding - shift_val4)];




figure;
set(gcf, 'color', 'w');
subplot(2,1,1);
plot(s1,'linewidth',3); box on; grid on; hold on;
plot(s1_1,'linewidth',2);
plot(s1_2,'linewidth',2);
plot(s1_3,'linewidth',2);
plot(s1_4,'linewidth',2);
xlabel('time (s)','Interpreter','Latex','FontSize',axis_font_size);
xticks = get(gca,'XTick')/Fs;
for k = 1:length(xticks)
    xticklabels{k} = num2str(xticks(k),3);
end
set(gca,'XTickLabels',xticklabels);
lgnd = legend('$s_{true1}$', '$s_{1,1}$','$s_{2,1}$','$s_{3,1}$','$s_{4,1}$','Interpreter','latex');
lgnd.FontSize= 16;
set(gcf,'Position',[1441 671 1218 649]);

subplot(2,1,2);
plot(s2,'linewidth',3); box on; grid on; hold on;
plot(s2_1,'linewidth',2);
plot(s2_2,'linewidth',2);
plot(s2_3,'linewidth',2);
plot(s2_4,'linewidth',2);
xlabel('time (s)','Interpreter','Latex','FontSize',axis_font_size);
xticks = get(gca,'XTick')/Fs;
for k = 1:length(xticks)
    xticklabels{k} = num2str(xticks(k),3);
end
set(gca,'XTickLabels',xticklabels);
lgnd = legend('$s_{true2}$', '$s_{1,2}$','$s_{2,2}$','$s_{3,2}$','$s_{4,2}$','Interpreter','latex');
lgnd.FontSize= 16;

if save_fig
    saveas(gcf, [figure_folder_path ex_name 'shifted_signals.png'])
end
%% Coherence calculations

% matrix1 = [s1_1 ;s2_1; s1_2; s2_2; s1_3; s2_3; s1_4; s2_4]';

matrix1 = [s1_1; s2_1 ; s2_2; s1_2;  s2_3; s1_3; s1_4; s2_4]';
%first source: 1,4,6,7
%second source:2,3,5,8

window = num_samples/5;
noverlap = window/2;
% f = 0:.25:5;
f = 0:5:100;
[cxy] = mscohere(s1,s2, window, noverlap, f, Fs)';
cxy(:,2) = f;

for i = 1:size(matrix1,2)
    for j = 1:size(matrix1,2)
        array_coh{i,j} = mscohere(matrix1(:,i), matrix1(:,j), window, noverlap, f, Fs);
    end
end

for k = 1:length(f)
    freq_index = find(f == f(k));
    for i = 1:size(array_coh,1)
        for j = 1:size(array_coh,2)
            array_coh_freq{k,1}(i,j) = array_coh{i,j}(freq_index);
        end
    end
end

%% Matrix Reordering

%Choose frequency for array coherence
freq = 10;
f_ind = find(f==freq);
coh_matrix = array_coh_freq{f_ind};


%Threshhold coherence matrix (40% of values)
B = coh_matrix;
threshold = prctile(B,40,'all');
B(B <= threshold) = 0;

idx_rcm = symrcm(B)';
idx_amd = symamd(B)';

coh_matrix_rcm = coh_matrix(idx_rcm,idx_rcm);
coh_matrix_amd = coh_matrix(idx_amd,idx_amd);

%% Plots for matrix reordering
figure; set(gcf, 'color', 'w');
imagesc(coh_matrix); box on; 
title('Coherence Matrix', 'fontsize', title_font_size);
colormap(copper); colorbar;

if save_fig
    saveas(gcf, [figure_folder_path ex_name 'coh_matrix.png'])
end

title_font_size = 20 ;
figure; 
subplot(1,2,1); set(gcf, 'color', 'w');
imagesc(coh_matrix_rcm); box on; 
sgtitle('Reordered Coherence Matrix', 'fontsize', title_font_size);
title('SRCM', 'fontsize', title_font_size);
colormap(copper); colorbar;

subplot(1,2,2); set(gcf, 'color', 'w');
imagesc(coh_matrix_amd); box on; 
title('SAMD', 'fontsize', title_font_size);
colormap(copper); colorbar;
set(gcf,'Position',[119 348 1244 449]);
if save_fig
    saveas(gcf, [figure_folder_path ex_name 'coh_matrix_reorder.png'])
end