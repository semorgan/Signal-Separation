%% Setup
clear; close all; clc; 

title_font_size = 14;
axis_font_size = 18;

fig_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Background Statistics/figures/';
save_fig = 1;
%% Visual Example of whitening data

rng(2); %set the seed for random number generator
mu = [0 0]; %mean
S = [1 .9; .9 3]; %covariance matrix
data = mvnrnd(mu, S, 1000)';

%Eigenvalue Decomposition
[E,D] = eig(S);
d_inv = diag(sqrt(1./diag(D)));

%Scale by D^(-1/2)
dataWhite = d_inv' * E * data;

%Plotting
figure('Position', [579 533 743 260]);
%Before Whitening
set(gcf, 'color', 'w');
subplot(1,2,1); hold on; grid on; box on;
title('Before Whitening', 'FontSize', title_font_size)
scatter(data(1,:),data(2,:),'filled');
xlabel('$\vec{x_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{x_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)

%After Whitening
subplot(1,2,2); hold on; grid on; box on;
title('After Whitening', 'FontSize', title_font_size)
scatter(dataWhite(1,:),dataWhite(2,:),'filled');
xlabel('$\vec{x_1}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)
ylabel('$\vec{x_2}$', 'Interpreter', 'latex', 'FontSize', axis_font_size)

%Save figure
if save_fig
   saveas(gcf, [fig_folder_path 'whitening.png']);
end
