%% Setup for running real data example
clear; close all; clc;


ex_name = 'magnetic_100_acoustic_280';

figure_folder_path = '/Users/semorgan/Documents/MATLAB/Github Folder/Signal-Separation/ICA/figures/writeup_figures/';
save_figures = 1;

data_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Test data/';
load([data_folder_path 'Acoustic280_Magnetic100.mat'])


axis_font_size = 16;
title_font_size = 14;




%% Acoustic and magnetic sources from lab testing

s1 = Dasdata_magnetic_100_acoustic_280.Strain(47,:);
s2 = Dasdata_magnetic_100_acoustic_280.Strain(42,:);

Fs =  1000/(Dasdata_magnetic_100_acoustic_280.Time(2)-Dasdata_magnetic_100_acoustic_280.Time(1));

%% Figures

figure; set(gcf, 'color', 'w');
subplot(2,1,1)
plot(s1);
xlim([0,length(s1)]);
ylim([min(s1),max(s1)])
xticks = get(gca,'XTick')/Fs;
for j = 1:length(xticks)
    xticklabels{j} = num2str(xticks(j),3)/1000;
end
set(gca,'XTickLabels',xticklabels);
xlabel('Time (sec)','Interpreter','Latex','FontSize',axis_font_size);
ylabel('Micro-strain ($\mu_\epsilon$)','Interpreter','Latex','FontSize',axis_font_size)
title('Channel 47 (acoustic and magnetic sensing fiber)','FontSize',title_font_size)
set(gcf,'color','w')


subplot(2,1,2)
plot(s2);
xlim([0,length(s2)]);
ylim([min(s2),max(s2)])
xticks = get(gca,'XTick')/Fs;
for j = 1:length(xticks)
    xticklabels{j} = num2str(xticks(j),3)/1000;
end
set(gca,'XTickLabels',xticklabels);
xlabel('Time (sec)','Interpreter','Latex','FontSize',axis_font_size);
ylabel('Micro-strain ($\mu_\epsilon$)','Interpreter','Latex','FontSize',axis_font_size)
title('Channel 42 (acoustic sensing fiber)','FontSize',title_font_size)
set(gcf,'color','w')


if save_figures
    saveas(gcf, [figure_folder_path, ex_name '_channel47_42.png'])
end

