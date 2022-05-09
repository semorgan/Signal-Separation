%% Setup for running real data example
clear; close all; clc;


ex_name = 'magnetic_100_acoustic_100';

figure_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/ICA/figures/';
save_fig = 0;

data_folder_path = '//Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Test data/';
load([data_folder_path 'Acoustic280.mat'])
load([data_folder_path 'Magnetic100.mat'])

axis_font_size = 16;
title_font_size = 14;


num_iter = 1000;

%% Acoustic and magnetic sources from lab testing

s1 = Dasdata_magnetic_100.Strain(47,:);
s2 = Dasdata_acoustic_280.Strain(47,:);

s1 = s1/norm(s1);
s2 = s2/norm(s2);

%make test cases the same length in time
s1(length(s2)+1:end) = [];


Fs =  1000/(Dasdata_acoustic_280.Time(2)-Dasdata_acoustic_280.Time(1));

true_mag1 = 100;
true_mag2 = 280;


S_orig = [s1; s2]; %each signal is a column



noise_level = 0;

%%
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
    % Mixtures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = rand(size(S_orig,1)); %mixing matrix
  
    X = A*S_orig;
    
    X = X + noise_level*randn(size(X));

    %% Preprocess the data and perform separation
    X_orig = X;

    [X,X_mean] = preprocess(X);

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
    [amp1_temp, locs1_temp] = findpeaks(S_fft(1,:), 'SortStr', 'descend');
    freq1_temp = freq(1,locs1_temp);
    locs1(1) = locs1_temp(min([find((mod(freq1_temp,true_mag1)<5) & (freq1_temp > true_mag1),1) find(abs(freq1_temp - true_mag1) < 5,1)]));
    locs1(2) = locs1_temp(min([find((mod(freq1_temp,true_mag2)<5) & (freq1_temp > true_mag2),1) find(abs(freq1_temp - true_mag2) < 5,1)]));
    amp1(1) = amp1_temp(locs1_temp == locs1(1));
    amp1(2) = amp1_temp(locs1_temp == locs1(2));
    newamps1 = [amp1(mod(locs1,true_mag1)<5) amp1(mod(locs1,true_mag2)<5)];
    
    
    [amp2_temp, locs2_temp] = findpeaks(S_fft(2,:), 'SortStr', 'descend');
    freq2_temp = freq(2,locs2_temp);
    locs2(1) = locs2_temp(min([find((mod(freq2_temp,true_mag1)<5) & (freq2_temp > true_mag1),1) find(abs(freq2_temp - true_mag1) < 5,1)]));
    locs2(2) = locs2_temp(min([find((mod(freq2_temp,true_mag2)<5) & (freq2_temp > true_mag2),1) find(abs(freq2_temp - true_mag2) < 5,1)]));
    amp2(1) = amp2_temp(locs2_temp == locs2(1));
    amp2(2) = amp2_temp(locs2_temp == locs2(2));
    newamps2 = [amp2(mod(locs2,true_mag1)<50) amp2(mod(locs2,true_mag2)<50)];
    
    
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
    %For the case when freq1 has 100 hz source 1st
    if (mod(freq1(1),true_mag1) < 5 && freq1(1) > true_mag1) || (abs(freq1(1) - true_mag1) < 5)
        if (mod(freq2(1),true_mag1) < 5 && freq2(1) > true_mag1) || (abs(freq2(1) - true_mag1) < 5)
            freq2(2) = temp_freq2(1);
            freq2(1) = temp_freq2(2);
            amp2(2) = temp_amp2(1);
            amp2(1) = temp_amp2(2);
        end
    %For the case when freq2 has 100 hz source 1st
    else
        if (mod(freq2(1),true_mag2) < 5 && freq2(1) > true_mag2) || (abs(freq2(1) - true_mag2) < 5)
            freq2(2) = temp_freq(1);
            freq2(1) = temp_freq(2);
            amp2(2) = temp_amp(1);
            amp2(1) = temp_amp(2);
        end
        
    end

    %Calculate relative magnitude between peak amplitudes
    mag1 = amp1(2)/amp1(1);
    mag2 = amp2(2)/amp2(1);


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
%% Run simulation_plots.m

if num_iter > 1
   
    simulation_plots
end




