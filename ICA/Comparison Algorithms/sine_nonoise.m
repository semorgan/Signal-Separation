%% Setup
clear;  clc;
close all;
format long e


ex_name = 'sine_100_40_nonoise';

figure_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/ICA/figures/';


save_fig = 0;
axis_font_size = 16;
title_font_size = 14;

% num_iter determine number of times mixing and separation is done
% set num_iter to 1 to do ICA one time for kurtosis and negentopy 
num_iter = 1000;



%% Example with sine waves at different frequencies
num_samples = 2000;
seconds = 3;
Fs = num_samples/seconds;

true_mag1 = 100;
true_mag2 = 40;

t = linspace(0,seconds,num_samples);
s1 = sin(2*pi*true_mag1*t);
s2 = sin(2*pi*true_mag2*t);

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
    
    X = X + noise_level*(2*rand(size(X)) - 1);
    
    %% Plot initial  data (true and mixed)
    axis_font_size = axis_font_size + 7;
    if num_iter == 1 && j == 1
        figure;
        
        set(gcf, 'color', 'w');
        set(gcf,'Position',[686 549 1222 633]);
        sgtitle('Source and Mixed Synthetic Signals','FontSize',title_font_size+7);
        subplot(4,6,1:5); 
        strue1 = plot(S_orig(1,:),'color',[0 0.4470 0.7410]); 
        ylabel('$s_{true1}$','Interpreter','Latex','FontSize',axis_font_size);
        xticks = get(gca,'XTick')/Fs;
        for k = 1:length(xticks)
            xticklabels{k} = num2str(xticks(k),3);
        end
        set(gca,'XTickLabels',[]);
        subplot(4,6,7:11); 
        strue2 = plot(S_orig(2,:),'color',[0.6350 0.0780 0.1840]);
        ylabel('$s_{true2}$','Interpreter','Latex','FontSize',axis_font_size);
        set(gca,'XTickLabels',[]);
        subplot(4,6,13:17); 
        mixed1 = plot(X(1,:),'color',[0.4940 0.1840 0.5560]);
        ylabel('$x_1$','Interpreter','Latex','FontSize',axis_font_size);
        set(gca,'XTickLabels',[]);
        subplot(4,6,19:23); 
        mixed2 = plot(X(2,:),'color',[0.4940 0.1840 0.5560]);
         ylabel('$x_2$','Interpreter','Latex','FontSize',axis_font_size);
        set(gca,'XTickLabels',xticklabels);
        xlabel('Time (sec)','FontSize',axis_font_size);
        lgnd = legend([strue1,strue2,mixed1],{'Source signal 1 (100 Hz)','Source signal 2 (40 Hz)','Mixed signals'});
        lgnd.FontSize= 16;
        lgnd.Position = [7.943535188216038e-01 4.541864139020536e-01 1.759410801963992e-01 1.216429699842021e-01]; 
        
        if save_fig
            saveas(gcf, [figure_folder_path ex_name '_source_mixed_signals.png'])
        end
    end
    axis_font_size = axis_font_size - 7;
    
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
    
    %If you are doing ICA one time and want to see the optimization plots
    if num_iter == 1 
        if j == 1
            kurt_plot = 1; %displays plots for kurtosis objective function
        else 
            kurt_plot = 0; %displays plots for negentropy objective function
        end
        optimization_plots
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
    %Update fields for each algorithm run
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
end
%% If doing separation more than 1 time, create bar graphs and tables
if num_iter > 1
    simulation_plots
end