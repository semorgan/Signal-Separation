%% Setup
close all; clear; clc


ex_name = 'matrix_reordering_example_structure';

figure_folder_path = '/Users/semorgan/Documents/MATLAB/Signal-Separation-Github/Array Coherence/figures/';


save_fig = 1;

axis_font_size = 16;
title_font_size = 14;

%% Matrix Reordering
B = bucky;

% Reverse Cuthill McKee algorithm
r = symrcm(B);
R = B(r,r);


% Symmetric Approximate Degree Permutation
p = symamd(B);
S = B(p,p);

% Plotting
figure;
set(gcf, 'color', 'w', 'Position',[441 581 578 216]);
sgtitle('Sparsity Structure','FontSize', title_font_size);
subplot(1,3,1),spy(B), title('Original Matrix','FontSize', title_font_size);box on;
subplot(1,3,2), spy(R,4), title('SRCM','FontSize', title_font_size); box on;
subplot(1,3,3), spy(S,4), title('SAMD','FontSize', title_font_size);box on;



if save_fig
    saveas(gcf, [figure_folder_path ex_name '.png'])
end



