%% Setup for running real data example
clear; close all; clc;

ex_name = 'coh_real_data_ac_guitar';
figure_folder_path = '/Users/semorgan/Documents/MATLAB/Github Folder/Signal-Separation/ICA/figures/writeup_figures/';
save_fig = 0;

axis_font_size = 16;
title_font_size = 14;

data_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Test data/';


%% Getting and formatting data

load([data_folder_path 'bass_guitar.mat'])
load([data_folder_path 'AC.mat'])

Fs = 62500;

num_samples = min(length(Dasdata_guitar.Strain), length(Dasdata_AC.Strain));

guitar = Dasdata_guitar.Strain(:,1:num_samples);
AC = Dasdata_AC.Strain(:,1:num_samples);


%Get 4 channels from each data set
s1 = guitar([20 25 30 35],:);
s2 = AC([20 25 30 35],:);


%% Coherence

matrix1 = [s1(1,:) ;s2(1,:); s1(2,:);s2(2,:);...
    s1(3,:); s2(3,:); s1(4,:); s2(4,:)]';

window = num_samples/5;
noverlap = window/2;
f = 0:5:100;


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

%% Reorder
%Choose frequency for array coherence
freq = 80;
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


%% Plots of coherence matrix


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


