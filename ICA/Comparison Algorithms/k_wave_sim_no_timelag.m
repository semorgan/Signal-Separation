%% Setup for running kwave simulation
clear; close all; clc;

ex_name = 'kwave_sine_60_25';
figure_folder_path = '/Users/semorgan/Documents/MATLAB/Github Folder/Signal-Separation/ICA/figures/writeup_figures/';

save_fig = 0;
title_font_size = 14;
axis_font_size = 16;

% num_iter determine number of times mixing and separation is done
% set num_iter to 1 to do ICA one time for kurtosis and negentopy 
num_iter = 50;

%% k-Wave setup
% create the computational grid
Nx = 128;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 0.25;        % grid point spacing in the x direction [m]
dy = 0.25;        % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);


%Define sources and sensors
src.locx = [Nx-Nx/4; Nx-Nx/4 ];
src.locy = [floor(Ny/3); ceil(2*Ny/3)];
src.freq = [60; 25];
src.mag = [1; 1];

sens.locx = [Nx/4; Nx/4];
sens.locy = [floor(Ny/3);ceil(2*Ny/3)];

% define the properties of the propagation medium
medium.sound_speed = 344; %1500;  % [m/s]
medium.alpha_coeff = 0.9; %0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.1; %1.5;

% create time array
t_end = 1;       % [s]
kgrid.makeTime(medium.sound_speed,[],t_end);

Fs = 1 / kgrid.dt;


%Define sources
source.p_mask = zeros(Nx, Ny);
for i = 1:length(src.locx)
    source.p_mask(src.locx(i), src.locy(i)) = 1;
    source.p(i,:) = src.mag(i) * sin(2 * pi * src.freq(i) * kgrid.t_array);
    % filter the source to remove high frequencies not supported by the grid
    source.p(i,:) = filterTimeSeries(kgrid, medium, source.p(i,:));
end

%Define sensor points
sensor.mask = zeros(Nx, Ny);
for i = 1:length(sens.locx)
    sensor.mask(sens.locx(i), sens.locy(i)) = 1;
end

%Plot source signals
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));
for i = 1:size(source.p,1)
    subplot(size(source.p,1),1,i)
    plot(kgrid.t_array * scale, source.p(i,:), 'k-');
    xlabel(['Time (' prefix 's)'],'FontSize',axis_font_size);
    ylabel('Signal Amplitude','FontSize',axis_font_size);
    axis tight;
    title('Input Pressure Signal','FontSize',title_font_size);
end
S_orig = source.p;

set(gcf, 'color', 'w'); 
if save_fig
    saveas(gcf, [figure_folder_path ex_name '_input_signals.png'])
end


% define the acoustic parameters to record
sensor.record = {'p', 'p_final'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);

% plot the final wave-field
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, ...
    sensor_data.p_final + source.p_mask + sensor.mask, [-1, 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
set(gcf, 'color', 'w'); 
if save_fig
    saveas(gcf, [figure_folder_path ex_name '_pressure_field.png'])
end

% plot the simulated sensor data
figure;
for i = 1:size(sensor_data.p,1)
    subplot(size(sensor_data.p,1), 1, i);
    plot(kgrid.t_array * scale, sensor_data.p(i,:), 'r-');
    xlabel(['Time (' prefix 's)'],'FontSize',axis_font_size);
    ylabel('Signal Amplitude','FontSize',axis_font_size);
    axis tight;
    title('Sensor Pressure Signal','FontSize',title_font_size);
end
set(gcf, 'color', 'w'); 
if save_fig
    saveas(gcf, [figure_folder_path ex_name '_resulting_signals.png'])
end

%% Perform separation 1000 times
for j = 1:2
    
    %Variables to add to on each loop through
    w1_initial = zeros(num_iter, 2);
    w2_initial = zeros(num_iter, 2);
    peak_freq1 = zeros(num_iter, 2);
    peak_freq2 = zeros(num_iter, 2);
    peak_amp1 = zeros(num_iter, 2);
    peak_amp2= zeros(num_iter, 2);
    rel_mag1 = zeros(num_iter, 1);
    rel_mag2 = zeros(num_iter, 1);
    
    for iter = 1:num_iter
    %% Preprocess the data and perform separation
    X_orig = sensor_data.p;

    [X,X_mean] = preprocess(X_orig);

    if j == 1
        %Fixed point iteration kurtosis algorithm
        [S,W, w1 ,w2] = FP_kurt(X,X_mean);
    else
        %Fast ICA Negentropy
        [S,W, w1 ,w2] = FP_negen(X,X_mean);
    end
    %% Run plots

    [freq, S_fft] = calc_fft(S, Fs);

    %Find amplitudes and frequencies for peaks after separation
    [amp1, locs1] = findpeaks(S_fft(1,:),'NPeaks',2, 'SortStr', 'descend');
    [amp2, locs2] = findpeaks(S_fft(2,:),'NPeaks',2, 'SortStr', 'descend');
    freq1 = freq(1,locs1);
    freq2 = freq(1,locs2);

    if length(freq1) == 1
        freq1(2) = 0;
        amp1(2) = 0;
    end
    if length(freq2) == 1
        freq2(2) = 0;
        amp2(2) = 0;
    end

    %Reorder frequency and amplitude based on peak amplitudes
    temp_amp2 = amp2;
    temp_freq2 = freq2;
    if freq1(1) <= freq2(1)
        freq2 = freq1;
        freq1 = temp_freq2;
        amp2 = amp1;
        amp1 = temp_amp2;
    end

    %Calculate relative magnitude between peak amplitudes
    mag1 = amp1(2)/amp1(1);
    mag2 = amp2(2)/amp1(1);


    %% Updating variables through each loop
    w1_initial(iter, :) = w1(:,1)';
    w2_initial(iter, :) = w2(:,1)'; 
    peak_freq1(iter, :) = freq1;
    peak_freq2(iter, :) = freq2;
    peak_amp1(iter, :) = amp1;
    peak_amp2(iter, :) = amp2;
    rel_mag1(iter, :) = mag1;
    rel_mag2(iter, :) = mag2;
    
    end
    %Update fields for each algorithm ran
    alg.w1_initial =  w1_initial;
    alg.w2_initial =  w2_initial;
    alg.peak_freq1 =  peak_freq1;
    alg.peak_freq2 =  peak_freq2;
    alg.peak_amp1 =  peak_amp1;
    alg.peak_amp2 =  peak_amp2;
    alg.rel_mag1 =  rel_mag1;
    alg.rel_mag2 =  rel_mag2;
    
    if j == 1
        FPkurt = alg;
        kurt_plot = true;
        negen_plot = false;
    else
        FPnegen = alg;
        negen_plot = true;
        kurt_plot =false;
    end
    
    clear alg;
    
   if num_iter == 1
       optimization_plots 
   end 
end

%%

if num_iter > 1
    
    true_mag1 = 60;
    true_mag2 = 25;
   
    simulation_plots 
end



