%% Setup for running real data example
clear; close all; clc;

figure_folder_path = '/Users/semorgan/Documents/MATLAB/Github Folder/Signal-Separation/ICA/figures/writeup_figures/';
% save_figures = 1;


%% Reading in Data

% Add files to current path
addpath('/Users/semorgan/Documents/MATLAB/DAS Data')  %Folder where all DAS Data is located
addpath('/Users/semorgan/Documents/MATLAB/Github Folder/Signal-Separation/DataProcessing/File Readers') %Folder with ReadDasLog versions

% If you want to manually select files write (won't work for saving figures)
% folder_path = []; file_start = [];

% Otherwise define path for folder containing data
Dasdata_folder = '/Users/semorgan/Documents/MATLAB/DAS Data/2022_03_18/';
testName = 'Acoustic_280_Magnetic_100'; %Use if multiple events stored in one folder


% file_start = 'Acoustic 280hz 100hz';
% plot_title = 'Acoustic (280 hz) Magnetic (100 hz)';


file_start = '100hz Mag'; 
plot_title = 'Magnetic (100 hz)';


% file_start = 'Acoustic 280 hz2';
% plot_title = 'Acoustic (280 hz)';

% if folder has single event put file_start = [];
save_plot = 1; %saves plots as png files if true
channels = [47]; %define specific channels you want plots of 

folder_path = [Dasdata_folder testName '/'];
figure_save_name = [testName '_' file_start];
Dasdata = ReadDasLogV5(folder_path, file_start);
% Fs = 62500; %sampling frequency

Fs =  1000/(Dasdata.Time(2)-Dasdata.Time(1));

initial_analysis

%Channels 47 (magnetic fiber) and 42 (sentek fiber)
% save('Acoustic280_Magnetic100','Dasdata','Fs')


Dasdata_magnetic_100 = Dasdata;
save('Magnetic100','Dasdata_magnetic_100')

% Dasdata_acoustic_280 = Dasdata;
% save('Acoustic280','Dasdata_acoustic_280')

% figure(1);
% ylim([41.5 47.5])
% if save_plot
%     saveas(gcf, [figure_folder_path, figure_save_name, '_Dasarray_channels.png'])
% end







