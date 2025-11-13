% clear
% clc

load('Boresight_MonteCarlo_100.mat')
lw = 3;
fs = 20;
fname = 'Times New Roman';
mss = 14;
afFigurePosition = [5 5 22 15]; 
psize = [10 10];
axesFontSize = fs;
legendFontSize = fs;
M_vec = [7];
plotvec = [];
figure
set(gcf, 'Units', 'centimeters'); 
set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
plotvec(end + 1) = plot(boresight_vec,Boresight_FD_Vec_SDP.*1e6, '-', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec,Boresight_FD_Vec_MRT.*1e6, '-.', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec,Boresight_FD_Vec_SDP.*1e6, '-s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_Line_Vec_SDP.*1e6, '-o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_Ploy_Vec_SDP.*1e6, '-*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_RS_Vec_SDP*1e6, '-d', 'color', "#77AC30", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_SGP_Vec_SDP*1e6, '-p', 'color', "#EDB120", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_SCA_Vec_SDP*1e6, '-h', 'color', "#7E2F8E", 'LineWidth', lw, 'MarkerSize',mss);
hold on
% Start main plot
plotvec(end + 1) = plot(boresight_vec,Boresight_FD_Vec_SDP.*1e6, '-s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_Line_Vec_SDP.*1e6, '-o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_Ploy_Vec_SDP.*1e6, '-*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Boresight_RS_Vec_SDP*1e6, '-d', 'color', "#77AC30", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec,Boresight_FD_Vec_MRT.*1e6, '-.s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Result_Line_Vec_MRT_overL.*1e6, '-.o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Result_Poly_Vec_MRT_overL.*1e6, '-.*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Result_RS_Vec_MRT_overL*1e6, '-.d', 'color', "#77AC30", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Result_SGP_Vec_MRT_overL.*1e6, '-.p', 'color', "#EDB120", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(boresight_vec ,Result_SCA_Vec_MRT_overL*1e6, '-.h', 'color', "#7E2F8E", 'LineWidth', lw, 'MarkerSize',mss);
hold on
box on
grid on
ylabel('Average $\min_i P_i^{rx}$ ($\mu$W)','Interpreter', 'latex','fontsize',axesFontSize)
xlabel('Radio stripe length (m)','Interpreter', 'latex','fontsize',axesFontSize)
xticks([2 4 8 16])
xticklabels({'2', '4', '8', '16'})
xlim([1.8 16.2])
set(gca,'FontSize',fs)
fontname(gcf, fname)
hl1 = legend([plotvec(1) plotvec(2)], 'SDP-based precoder', 'MRT-based precoder', 'Location', 'northwest');
ah2 = axes('position',get(gca,'position'),'visible','off');
hl2 = legend(ah2, [plotvec(3) plotvec(4) plotvec(5) plotvec(6) plotvec(7) plotvec(8)], ...
    'Center-FD', 'SGP-Line', 'SGP-Polygon', 'Center-Rectangle', 'SGP', 'SCA', 'Location', 'northwest');
set([hl1 hl2],'interpreter','latex','FontSize',legendFontSize);

% saveas(gcf, strcat('SDP_F', string(ftest), '.fig'))
% saveas(gcf, strcat('SDP_F', string(ftest), '.eps'))

