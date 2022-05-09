%% Setup
clear; close all; clc;


ex_name = 'coh_square_hat_noshift';
figure_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Array Coherence/figures/';

save_fig = 0;

axis_font_size = 16;
title_font_size = 14;

%% Square and hat functions

a = 0;
b = 1;
A1 = 1;
A2 = 1;


Fs = 50000;
num_samples = Fs*(b-a);
t = linspace(a,b,Fs*(b-a));


padding = Fs/2;


s1 = A1*ones(size(t));
s2_1 = (2*A2*t(1:length(t)/2)/(b-a))+(2*A2*a)/(a-b);
s2_2 =  (2*A2*t(length(t)/2+1:end)/(a-b))+(2*A2*b)/(b-a);
s1 = [zeros(1,padding) s1 zeros(1,padding)];
s2 = [zeros(1,padding) s2_1 s2_2 zeros(1,padding)];


figure;
set(gcf, 'color', 'w');
% subplot(2,1,1); 
plot(s1,'color',[0 0.4470 0.7410],'linewidth',2); box on; grid on; hold on;
xlabel('time (s)','Interpreter','Latex','FontSize',axis_font_size);
xticks = get(gca,'XTick')/Fs;
for k = 1:length(xticks)
    xticklabels{k} = num2str(xticks(k),3);
end
set(gca,'XTickLabels',xticklabels);
plot(s2,'color',[0.8500 0.3250 0.0980],'linewidth',2);
lgnd = legend('$s_1$','$s_2$','Interpreter','latex');
lgnd.FontSize= 16;
set(gcf,'position',[1441 1024 645 296])

if save_fig
    saveas(gcf, [figure_folder_path ex_name '.png'])
end
%% Coherence calculations

matrix1 = [s1 ;s2; s2; s1; s2; s1; s1; s2]';

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


