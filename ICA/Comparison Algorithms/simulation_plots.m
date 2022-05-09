%% Run all plots for simulation (of 1000 algorithm runs)


figure; set(gcf, 'color', 'w');
sgtitle('Frequency of peak amplitude (estimated sources)','FontSize',title_font_size)
subplot(211); hold on; grid on; box on;
ylabel('Source 1: frequency (Hz)','FontSize', axis_font_size)
xlabel('iterations','FontSize', axis_font_size)
plot(1:num_iter,FPkurt.peak_freq1(:,1),'-o');
plot(1:num_iter,FPnegen.peak_freq1(:,1),'-o');
legend('Fixed point Kurtosis', 'Fixed point Negentropy')
subplot(212); hold on; grid on; box on;
ylabel('Source 2: frequency (Hz)','FontSize', axis_font_size)
xlabel('run number','FontSize', axis_font_size)
plot(1:num_iter,FPkurt.peak_freq2(:,1),'-o');
plot(1:num_iter,FPnegen.peak_freq2(:,1),'-o');
legend('Fixed point Kurtosis', 'Fixed point Negentropy')

if save_fig
    saveas(gcf, [figure_folder_path ex_name '_peak_amp_freq.png'])
end

 %Relative magnitude of peak frequency amplitudes
figure('Position', [689 419 948 536], 'color', 'w'); 
sgtitle('Relative Magnitude of Peak Frequency Amplitudes','FontSize',title_font_size)

subplot(221); box on; grid on; hold on;
histogram(FPkurt.rel_mag1,'BinEdges',(0:.05:1),'FaceColor',[0 0.4470 0.7410],'Normalization','probability')
ylabel('Source 1 (% algorithm runs)','FontSize', axis_font_size)
title('Fixed Point Kurtosis','FontSize',title_font_size)
ylim([0 1])
yticklabels(yticks*100)
subplot(222); box on; grid on; hold on;
histogram(FPnegen.rel_mag1,'BinEdges',(0:.05:1),'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability')
title('Fixed Point Negentropy','FontSize',title_font_size)
ylim([0 1])
yticklabels(yticks*100)

 subplot(223); box on; grid on; hold on;
histogram(FPkurt.rel_mag2,'BinEdges',(0:.05:1),'FaceColor',[0 0.4470 0.7410],'Normalization','probability')
ylabel('Source 2 (% algorithm runs)','FontSize', axis_font_size)
%     xlabel('$\big|\hat{s_1}(\omega_2^{(1)})/\hat{s_1}(\omega_1^{(1)})\big|$','Interpreter', 'Latex')
xlabel('Ratio of 2nd highest to highest peak amplitude','FontSize', axis_font_size);
ylim([0 1])
yticklabels(yticks*100)
subplot(224); box on; grid on; hold on;
histogram(FPnegen.rel_mag2,'BinEdges',(0:.05:1),'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability')
%     xlabel('$\big|\hat{s_2}(\omega_2^{(2)})/\hat{s_2}(\omega_1^{(2)})\big|$','Interpreter', 'Latex')
xlabel('Ratio of 2nd highest to highest peak amplitude','FontSize', axis_font_size);
ylim([0 1])
yticklabels(yticks*100)


if save_fig
    saveas(gcf, [figure_folder_path ex_name '_rel_mag_peak_freq_all.png'])
end

%% Calculate ratio totals for tables in appendix

% Kurtosis 1st source
FPkurt.ratio1(1) = length(find(FPkurt.rel_mag1<=0.1));
FPkurt.ratio1(2) = length(find(FPkurt.rel_mag1 <= 0.2 & FPkurt.rel_mag1>0.1));
FPkurt.ratio1(3) = length(find(FPkurt.rel_mag1 <= 0.3 & FPkurt.rel_mag1>0.2));
FPkurt.ratio1(4) = length(find(FPkurt.rel_mag1 <= 0.4 & FPkurt.rel_mag1>0.3));
FPkurt.ratio1(5) = length(find(FPkurt.rel_mag1 <= 0.5 & FPkurt.rel_mag1>0.4));
FPkurt.ratio1(6) = length(find(FPkurt.rel_mag1 <= 0.6 & FPkurt.rel_mag1>0.5));
FPkurt.ratio1(7) = length(find(FPkurt.rel_mag1 <= 0.7 & FPkurt.rel_mag1>0.6));
FPkurt.ratio1(8) = length(find(FPkurt.rel_mag1 <= 0.8 & FPkurt.rel_mag1>0.7));
FPkurt.ratio1(9) = length(find(FPkurt.rel_mag1 <= 0.9 & FPkurt.rel_mag1>0.8));
FPkurt.ratio1(10) = length(find(FPkurt.rel_mag1 <= 1 & FPkurt.rel_mag1>0.9));

% Kurtosis 2nd source
FPkurt.ratio2(1) = length(find(FPkurt.rel_mag2<=0.1));
FPkurt.ratio2(2) = length(find(FPkurt.rel_mag2 <= 0.2 & FPkurt.rel_mag2>0.1));
FPkurt.ratio2(3) = length(find(FPkurt.rel_mag2 <= 0.3 & FPkurt.rel_mag2>0.2));
FPkurt.ratio2(4) = length(find(FPkurt.rel_mag2 <= 0.4 & FPkurt.rel_mag2>0.3));
FPkurt.ratio2(5) = length(find(FPkurt.rel_mag2 <= 0.5 & FPkurt.rel_mag2>0.4));
FPkurt.ratio2(6) = length(find(FPkurt.rel_mag2 <= 0.6 & FPkurt.rel_mag2>0.5));
FPkurt.ratio2(7) = length(find(FPkurt.rel_mag2 <= 0.7 & FPkurt.rel_mag2>0.6));
FPkurt.ratio2(8) = length(find(FPkurt.rel_mag2 <= 0.8 & FPkurt.rel_mag2>0.7));
FPkurt.ratio2(9) = length(find(FPkurt.rel_mag2 <= 0.9 & FPkurt.rel_mag2>0.8));
FPkurt.ratio2(10) = length(find(FPkurt.rel_mag2 <= 1 & FPkurt.rel_mag2>0.9));

%Negentropy 1st source
FPnegen.ratio1(1) = length(find(FPnegen.rel_mag1<=0.1));
FPnegen.ratio1(2) = length(find(FPnegen.rel_mag1 <= 0.2 & FPnegen.rel_mag1>0.1));
FPnegen.ratio1(3) = length(find(FPnegen.rel_mag1 <= 0.3 & FPnegen.rel_mag1>0.2));
FPnegen.ratio1(4) = length(find(FPnegen.rel_mag1 <= 0.4 & FPnegen.rel_mag1>0.3));
FPnegen.ratio1(5) = length(find(FPnegen.rel_mag1 <= 0.5 & FPnegen.rel_mag1>0.4));
FPnegen.ratio1(6) = length(find(FPnegen.rel_mag1 <= 0.6 & FPnegen.rel_mag1>0.5));
FPnegen.ratio1(7) = length(find(FPnegen.rel_mag1 <= 0.7 & FPnegen.rel_mag1>0.6));
FPnegen.ratio1(8) = length(find(FPnegen.rel_mag1 <= 0.8 & FPnegen.rel_mag1>0.7));
FPnegen.ratio1(9) = length(find(FPnegen.rel_mag1 <= 0.9 & FPnegen.rel_mag1>0.8));
FPnegen.ratio1(10) = length(find(FPnegen.rel_mag1 <= 1 & FPnegen.rel_mag1>0.9));

%Negentropy 2nd source
FPnegen.ratio2(1) = length(find(FPnegen.rel_mag2<=0.1));
FPnegen.ratio2(2) = length(find(FPnegen.rel_mag2 <= 0.2 & FPnegen.rel_mag2>0.1));
FPnegen.ratio2(3) = length(find(FPnegen.rel_mag2 <= 0.3 & FPnegen.rel_mag2>0.2));
FPnegen.ratio2(4) = length(find(FPnegen.rel_mag2 <= 0.4 & FPnegen.rel_mag2>0.3));
FPnegen.ratio2(5) = length(find(FPnegen.rel_mag2 <= 0.5 & FPnegen.rel_mag2>0.4));
FPnegen.ratio2(6) = length(find(FPnegen.rel_mag2 <= 0.6 & FPnegen.rel_mag2>0.5));
FPnegen.ratio2(7) = length(find(FPnegen.rel_mag2 <= 0.7 & FPnegen.rel_mag2>0.6));
FPnegen.ratio2(8) = length(find(FPnegen.rel_mag2 <= 0.8 & FPnegen.rel_mag2>0.7));
FPnegen.ratio2(9) = length(find(FPnegen.rel_mag2 <= 0.9 & FPnegen.rel_mag2>0.8));
FPnegen.ratio2(10) = length(find(FPnegen.rel_mag2 <= 1 & FPnegen.rel_mag2>0.9));

%Making table of these values
Ratio = {'0-0.1';'01.-0.2';'0.2-0.3';'0.3-0.4';'0.4-0.5'; '0.5-0.6'; '0.6-0.7';...
    '0.7-0.8'; '0.8-0.9'; '0.9-1'};
Kurt_src1 = FPkurt.ratio1';
Kurt_src2 = FPkurt.ratio2';
Negen_src1 = FPnegen.ratio1';
Negen_src2 = FPnegen.ratio2';
T = table(Kurt_src1,Kurt_src2, Negen_src1, Negen_src2,'RowNames',Ratio);

%Turning table into figure for saving purposes
figure;
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);


if save_fig
    saveas(gcf, [figure_folder_path ex_name '_rel_mag_peak_freq_table.png'])
end
