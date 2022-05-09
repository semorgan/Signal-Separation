%Script creates plots for analysis of optimization


if save_fig
    if kurt_plot 
        figure_save_name = [figure_folder_path ex_name '_FPkurt' ];
    else
        figure_save_name = [figure_folder_path ex_name '_FPnegen' ];
    end
end

title_font_size = 14;
axis_font_size = 16;


%% Plot mixed components

figure;
set(gcf, 'color', 'w');
sgtitle('Observed Mixed Components','FontSize',title_font_size)
subplot(2,1,1);  hold on; grid on; box on;
plot(X_orig(1,:));
title('$\vec{s_1}$','Interpreter','latex','FontSize',title_font_size)
subplot(2,1,2); hold on; grid on; box on;
plot(X_orig(2,:));
title('$\vec{s_2}$','Interpreter','latex','FontSize',title_font_size)



%% Plotting Normalized Independent components after separation
figure;
set(gcf, 'color', 'w');
sgtitle('Normalized Independent Components (estimated)','FontSize',title_font_size)
subplot(2,1,1); hold on; grid on; box on;
plot(S(1,:)/max(abs(S(1,:))));
subplot(2,1,2); hold on; grid on; box on;
plot(S(2,:)/max(abs(S(2,:))));


%Save figure
if save_fig
    saveas(gcf, [figure_save_name '_indep_comp_sep.png']);
end



%% Plot kurtosis as function of angle (optimization landscape)

%Find kurtosis values for angles from 0 to 2*pi
angle = linspace(0,2*pi,1000);
kurt = zeros(size(angle));
negen = zeros(size(angle));
x = cos(angle);
y = sin(angle);
w_angle = [x; y];

G = @(x) log(cosh(x));

for i = 1:length(angle)
    kurt(i) = abs(kurtosis(w_angle(:,i)'*X)-3);
    negen(i) = mean(G(w_angle(:,i)'*X),2);
end
%Find angle of 1st row of W matrix
w1_angle = atan2(w1(2,1:end),w1(1,1:end));
for i = 1:length(w1_angle)
    if w1_angle(i) < 0 
        w1_angle(i) = w1_angle(i) + 2*pi;
    end
    [~,w1_index(i)] = min(abs(angle-w1_angle(i)));
end
%Find angle of 2nd row of W matrix
w2_angle = atan2(w2(2,1:end),w2(1,1:end));
for i = 1:length(w2_angle)
    if w2_angle < 0
        w2_angle = w2_angle + 2*pi;
    end
    [~,w2_index(i)] = min(abs(angle-w2_angle(i)));
end

figure; set(gcf, 'color', 'w');
set(gcf,'Position',[441 379 935 420])
subplot(1,4,1:3); grid on; box on; hold on;
if kurt_plot
    %Creating plot for kurtosis
    title('Kurtosis Optimization Landscape','FontSize',title_font_size)
    plot(angle,kurt)
    ylim([(min(kurt)-(max(kurt)-min(kurt))*.1) (max(kurt)+(max(kurt)-min(kurt))*.1)])
    xlabel('angle of $w$ (rad)', 'Interpreter','latex','FontSize',axis_font_size)
    ylabel('$|$Kurtosis$|$','FontSize',axis_font_size,'Interpreter','latex')
    colors = get(gca,'colororder');
    plot(angle(w1_index),kurt(w1_index),'-o','MarkerFaceColor',colors(2,:),'MarkerSize',8)
    plot(angle(w2_index),kurt(w2_index),'-o','MarkerFaceColor',colors(3,:),'MarkerSize',8)
    scatter(angle(w1_index(end)),kurt(w1_index(end)),400,'*','k')
    scatter(angle(w2_index(end)),kurt(w2_index(end)),400,'*','k')
    lgnd = legend('$kurt(\vec{w}^\top\mathbf{x}$)','$\vec{w_1}^\top$','$\vec{w_2}^\top$','interpreter','latex');
    lgnd.FontSize = 16;

else 
    %Creating plot for negentropy
    title('Negentropy Optimization Landscape','FontSize',title_font_size)
    plot(angle,negen)
    xlabel('angle of $w$ (rad)', 'Interpreter','latex','FontSize',axis_font_size)
    ylabel('Negentropy','FontSize',axis_font_size, 'Interpreter','latex')
    colors = get(gca,'colororder');
    plot(angle(w1_index),negen(w1_index),'-o','MarkerFaceColor',colors(2,:),'MarkerSize',8)
    plot(angle(w2_index),negen(w2_index),'-o','MarkerFaceColor',colors(3,:),'MarkerSize',8)
    scatter(angle(w1_index(end)),negen(w1_index(end)),400,'*','k')
    scatter(angle(w2_index(end)),negen(w2_index(end)),400,'*','k')
    lgnd = legend('$E[G(\vec{w}^\top\mathbf{x}$)]','$\vec{w_1}^\top$','$\vec{w_2}^\top$','interpreter','latex');
    lgnd.FontSize = 16;
    
end
%Plotting  unit circle with w iterations
subplot(1,4,4); grid on; box on; hold on;
title('Weight Vectors','FontSize',title_font_size)

% draw unit circle
rectangle('Position',[-1 -1 2 2],'Curvature',[1 1],'Linewidth',1.5)
axis equal

%Plotting iterations of finding W matrix
plot([],[]);
plot(w1(1,:),w1(2,:),'-o','Color',colors(2,:),'MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:),'MarkerSize',8);
plot(w2(1,:),w2(2,:),'-o','Color',colors(2,:),'MarkerFaceColor',colors(3,:),'MarkerEdgeColor',colors(3,:),'MarkerSize',8);
scatter(w1(1,end),w1(2,end),400,'*','k')
scatter(w2(1,end),w2(2,end),400,'*','k')
ylim([-1.5 1.5]); xlim([-1.5 1.5]);
lgnd = legend('$\vec{w_1}$','$\vec{w_2}$','Interpreter','Latex');
lgnd.FontSize = 16;
lgnd.Location = 'northoutside';


%Save figure
if save_fig
    saveas(gcf, [figure_save_name '_opt_landscape.png']);
end
    


%% FFT plot
figure; set(gcf, 'color', 'w'); 
fig_fft = get(gcf,'Number');
    
    for i = 1:size(S,1)
        [f, P1] = calc_fft(S(i,:), Fs);
     
        figure(fig_fft);
        subplot(size(S,1),1,i); grid on; box on; hold on;
%         semilogy(f,P1) %y axis in log scale
        plot(f,P1)
        ylabel('Amplitude','FontSize',axis_font_size)
        
        S_fft(i,:) = P1;
        freq(i,:) = f;
        
    end
    figure(fig_fft)
    sgtitle('Single-Sided Amplitude Spectrum (estimated)','FontSize',title_font_size)
    xlabel('f (Hz)','FontSize',axis_font_size)
    
    if save_fig
       saveas(gcf, [figure_save_name, '_fft.png'])
    end
    
%% FFT plot original data
figure; set(gcf, 'color', 'w'); 
fig_fft = get(gcf,'Number');
    
    for i = 1:size(S_orig,1)
        
        [f, tP1] = calc_fft(S_orig(i,:), Fs);
    
        figure(fig_fft);
        subplot(size(S_orig,1),1,i); grid on; box on; hold on;
        plot(f,P1)
        if i == 1
            ylabel('Amplitude ($s_1$)','FontSize',axis_font_size, 'Interpreter', 'latex');
        else
             ylabel('Amplitude ($s_2$)','FontSize',axis_font_size, 'Interpreter', 'latex')
        end
        xlim([0 600])
        S_fft(i,:) = P1;
        freq(i,:) = f;
        
    end
    figure(fig_fft)
    sgtitle('Single-Sided Amplitude Spectrum','FontSize',title_font_size)
    xlabel('f (Hz)','FontSize',axis_font_size)
    
    if save_fig
       saveas(gcf, [figure_save_name, '_fft_originaldata.png'])
    end
    
    
 %% FFT mixed data
figure; set(gcf, 'color', 'w'); 
fig_fft = get(gcf,'Number');
    
    for i = 1:size(X,1)
        
        [f, P1] = calc_fft(X(i,:), Fs);
        
        figure(fig_fft);
        subplot(size(S_orig,1),1,i); grid on; box on; hold on;
        plot(f,P1)
        if i == 1
            ylabel('Amplitude ($s_1$)','FontSize',axis_font_size, 'Interpreter', 'latex');
        else
             ylabel('Amplitude ($s_2$)','FontSize',axis_font_size, 'Interpreter', 'latex')
        end
        xlim([0 600])
        S_fft(i,:) = P1;
        freq(i,:) = f;
        
    end
    figure(fig_fft)
    sgtitle('Single-Sided Amplitude Spectrum (mixed data)','FontSize',title_font_size)
    xlabel('f (Hz)','FontSize',axis_font_size)
    
    if save_fig
       saveas(gcf, [figure_save_name, '_fft_mixeddata.png'])
    end
    

    
   
