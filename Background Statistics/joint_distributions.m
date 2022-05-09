%% Setup
clear; close all; clc;

title_font_size = 14;
axis_font_size = 18;

fig_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Background Statistics/figures/';
save_fig = 1;

%% Creating random variables with different distributions

%Non - gaussian uniform components (individual)
s1 = rand(1000,1)';
s2 = rand(1000,1)';
s = [s1; s2];

%Gaussian components (individual)
s1g = normrnd(0, 1, 1000, 1)';
s2g = normrnd(0, 1, 1000, 1)';
sg = [s1g; s2g];

%Orthogonal mixing matrix & mixed signals
A = 1/sqrt(2)*[1 1; 1 -1];
x = A*s;
xg = A*sg;

%% Create figures showing joint distributions (2 figures per plot)

% non-Gaussian independent components
figure('Position', [1542 888 713 271]);
set(gcf, 'color', 'w');
subplot(1,2,1); hold on; grid on; box on;
scatter(s1, s2, 'filled')
xlabel('$\vec{s_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{s_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
title('Independent Components', 'FontSize', title_font_size)

% non-Gaussian mixed components
subplot(1,2,2); hold on; grid on; box on;
scatter(x(1,:), x(2,:), 'filled')
xlabel('$\vec{x_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{x_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
title('Mixed Components', 'FontSize', title_font_size)
%Save figure
if save_fig
    saveas(gcf, [fig_folder_path 'distr_nongauss.png']);
end

% Gaussian independent components
figure('Position', [1542 888 713 271]);
set(gcf, 'color', 'w');
subplot(1,2,1); hold on; grid on; box on;
scatter(s1g, s2g, 'filled')
xlabel('$\vec{s_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{s_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
title('Independent Components', 'FontSize', title_font_size)

% Gaussian mixed components
subplot(1,2,2); hold on; grid on; box on;
scatter(xg(1,:), xg(2,:), 'filled')
xlabel('$\vec{x_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{x_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
title('Mixed Components', 'FontSize', title_font_size)
%Save fugure
if save_fig
    saveas(gcf, [fig_folder_path 'distr_gauss.png']);
end


%% Calculate kurtosis and negentropy for different distributions


N=10000;
figure; 
set(gcf, 'color', 'w');

%Square wave
T=64;
yS=double(2*mod(int16(2*(1:N)/T),2)-1); %square wave, mean=0, stdev=1, period=T
subplot(3,1,1);         
hS=histogram(yS,'Normalization','pdf','BinWidth',.01);
title('Square wave','FontSize', title_font_size);   %plot title
entropyS=-hS.BinWidth*sum(hS.Values(hS.Values>0).*log(hS.Values(hS.Values>0))/log(2));
JS=log(std(yS))/log(2)+2.0471-entropyS;
kurtS = kurtosis(yS) - 3;


%Uniform noise
yU=(rand(1,N)-.5)*sqrt(12);     %uniform random noise, mean=0, stdev=1
subplot(3,1,2);         
hU=histogram(yU,'Normalization','pdf');
title('Uniform Distribution','FontSize', title_font_size);
entropyU=-hU.BinWidth*sum(hU.Values(hU.Values>0).*log(hU.Values(hU.Values>0))/log(2));
JU=log(std(yU))/log(2)+2.0471-entropyU;
kurtU = kurtosis(yU) - 3;


%Gaussian noise
yG=randn(1,N);      %Gaussian random noise, mean=0, stdev=1
subplot(3,1,3);    
hG=histogram(yG,'Normalization','pdf');
title('Gaussian Distribution','FontSize', title_font_size);
entropyG=-hG.BinWidth*sum(hG.Values(hG.Values>0).*log(hG.Values(hG.Values>0))/log(2));
JG=log(std(yG))/log(2)+2.0471-entropyG;
kurtG = kurtosis(yG) - 3;


%Save figure
if save_fig
    saveas(gcf, [fig_folder_path 'distr_kurt_negen.png']);
end


