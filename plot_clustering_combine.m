clear
clc

load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_2.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_3.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_4.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_5.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_6.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_7.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_8.mat')
load('C:\Users\aazarbah23\Desktop\PhD\Simulation\BS_Deploy_final\WCNC\Cluster_MonteCarlo_1.mat')


Combine_Result_FD_Vec_MRT_CL = [Combine_Result_FD_Vec_MRT_CL1 min(Combine_Result_FD_Vec_MRT_CL2) ...
    min(Combine_Result_FD_Vec_MRT_CL3) min(Combine_Result_FD_Vec_MRT_CL4) ...
    min(Combine_Result_FD_Vec_MRT_CL5)  ...
    min(Combine_Result_FD_Vec_MRT_CL6) min(Combine_Result_FD_Vec_MRT_CL7) ...
    min(Combine_Result_FD_Vec_MRT_CL8)];
Combine_Result_FD_Vec_SDP_CL = [Combine_Result_FD_Vec_SDP_CL1 min(Combine_Result_FD_Vec_SDP_CL2) ...
    min(Combine_Result_FD_Vec_SDP_CL3) min(Combine_Result_FD_Vec_SDP_CL4) ...
    min(Combine_Result_FD_Vec_SDP_CL5)  ...
    min(Combine_Result_FD_Vec_SDP_CL6) min(Combine_Result_FD_Vec_SDP_CL7) ...
    min(Combine_Result_FD_Vec_SDP_CL8)];
Combine_Result_poly_Vec_SDP_CL = [Combine_Result_poly_Vec_SDP_CL1 min(Combine_Result_poly_Vec_SDP_CL2) ...
    min(Combine_Result_poly_Vec_SDP_CL3) min(Combine_Result_poly_Vec_SDP_CL4) ...
    min(Combine_Result_poly_Vec_SDP_CL5)  ...
    min(Combine_Result_poly_Vec_SDP_CL6) min(Combine_Result_poly_Vec_SDP_CL7) ...
    min(Combine_Result_poly_Vec_SDP_CL8)];
Combine_Result_poly_Vec_MRT_CL = [Combine_Result_poly_Vec_MRT_CL1 min(Combine_Result_poly_Vec_MRT_CL2) ...
    min(Combine_Result_poly_Vec_MRT_CL3) min(Combine_Result_poly_Vec_MRT_CL4) ...
    min(Combine_Result_poly_Vec_MRT_CL5)  ...
    min(Combine_Result_poly_Vec_MRT_CL6) min(Combine_Result_poly_Vec_MRT_CL7) ...
    min(Combine_Result_poly_Vec_MRT_CL8)];
Combine_Result_SCA_Vec_MRT_CL = [Combine_Result_SCA_Vec_MRT_CL1 min(Combine_Result_SCA_Vec_MRT_CL2) ...
    min(Combine_Result_SCA_Vec_MRT_CL3) min(Combine_Result_SCA_Vec_MRT_CL4) ...
    min(Combine_Result_SCA_Vec_MRT_CL5)  ...
    min(Combine_Result_SCA_Vec_MRT_CL6) min(Combine_Result_SCA_Vec_MRT_CL7) ...
    min(Combine_Result_SCA_Vec_MRT_CL8)];
Combine_Result_SCA_Vec_SDP_CL = [Combine_Result_SCA_Vec_SDP_CL1 min(Combine_Result_SCA_Vec_SDP_CL2) ...
    min(Combine_Result_SCA_Vec_SDP_CL3) min(Combine_Result_SCA_Vec_SDP_CL4) ...
    min(Combine_Result_SCA_Vec_SDP_CL5)  ...
    min(Combine_Result_SCA_Vec_SDP_CL6) min(Combine_Result_SCA_Vec_SDP_CL7) ...
    min(Combine_Result_SCA_Vec_SDP_CL8)];



lw = 2.5;
fs = 14;
fname = 'Times New Roman';
mss = 8;
afFigurePosition = [5 5 15 10]; 
psize = [10 10];
axesFontSize = fs;
legendFontSize = fs;
M_vec = [7];
plotvec = [];
figure
set(gcf, 'Units', 'centimeters'); 
set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
plotvec(end + 1) = plot(ET_vec,pow2db(Combine_Result_SCA_Vec_SDP_CL), '-', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec,pow2db(Combine_Result_SCA_Vec_MRT_CL), '--', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec,pow2db(Combine_Result_FD_Vec_MRT_CL), '--s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_poly_Vec_MRT_CL), '--o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_SCA_Vec_MRT_CL), '--*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
% Start main plot
plotvec(end + 1) = plot(ET_vec,pow2db(Combine_Result_FD_Vec_SDP_CL), '-s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_poly_Vec_SDP_CL), '-o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_SCA_Vec_SDP_CL), '-*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec,pow2db(Combine_Result_FD_Vec_MRT_CL), '--s', 'color', "#0072BD", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_poly_Vec_MRT_CL), '--o', 'color', "#A2142F", 'LineWidth', lw, 'MarkerSize',mss);
hold on
plotvec(end + 1) = plot(ET_vec ,pow2db(Combine_Result_SCA_Vec_MRT_CL), '--*', 'color', "k", 'LineWidth', lw, 'MarkerSize',mss);
hold on
box on
grid on
ylabel('Average $\min_i P_i^{rx}$ ($\mu$W)','Interpreter', 'latex','fontsize',axesFontSize)
xlabel('$K$','Interpreter', 'latex','fontsize',axesFontSize)
% xticks([0.5 1 1.5 2 2.5 3])
% xticklabels({'0.5', '1', '1.5', '2', '2.5', '3'})
set(gca,'FontSize',fs)
fontname(gcf, fname)
hl1 = legend([plotvec(1) plotvec(2)], 'SDP-based precoder', 'MRT-based precoder', 'Location', 'northwest');
ah2 = axes('position',get(gca,'position'),'visible','off');
hl2 = legend(ah2, [plotvec(3) plotvec(4) plotvec(5)], ...
    'Center-FD','SGP-Polygon', 'SCA', 'Location', 'northwest');
set([hl1 hl2],'interpreter','latex','FontSize',legendFontSize);
set([hl1 hl2],'color','none');
set([hl1 hl2],'edgecolor','none');