clear
clc

Fvec = 2.5e9:2.5e9:30e9;
Nvec = [9 16 25 36 49 64 81 100 121 144 169 196];

Lvec = zeros(length(Fvec), length(Nvec));
fran_vec_poly = zeros(length(Fvec), length(Nvec));
fres_vec_poly = zeros(length(Fvec), length(Nvec));
fran_vec_line = zeros(length(Fvec), length(Nvec));
fres_vec_line = zeros(length(Fvec), length(Nvec));
fran_vec_FD = zeros(length(Fvec), length(Nvec));
fres_vec_FD = zeros(length(Fvec), length(Nvec));

for ff = 4:4
    f0 = Fvec(ff);
    lambda = 3e8/f0;
    inter_dist = lambda/2;
    for nn = 1:length(Nvec)
        N = Nvec(nn);
        maxL = (N - 1)*inter_dist;
        Lvec(ff, nn) = maxL;
        D = sqrt(2)*(sqrt(N) - 1)*inter_dist;
        fres_vec_FD(ff, nn) = ((D^4)/(8*lambda))^(1/3);
        fran_vec_FD(ff, nn) = (2*(D^2))/lambda;
        D = maxL;
        fres_vec_line(ff, nn) = ((D^4)/(8*lambda))^(1/3);
        fran_vec_line(ff, nn) = (2*(D^2))/lambda;
        r0 = inter_dist/(2*sin(pi/N));
        x = zeros(N, 2);
        phii = (2*pi)/N;
        for i = 1:N
            x(i, :) = [r0*cos((i-1)*phii) r0*sin((i-1)*phii)];
        end
        dvec = zeros(N, N);
        for i = 1:N
            for j = 1:N
                dvec(i, j) = norm(x(i, :) - x(j, :));
            end
        end
        D = max(max(dvec));
        fres_vec_poly(ff, nn) = ((D^4)/(8*lambda))^(1/3);
        fran_vec_poly(ff, nn) = (2*(D^2))/lambda;
    end
end

lw = 3;
fs = 16;
fname = 'Times New Roman';
mss = 12;
afFigurePosition = [5 5 25 12]; 
psize = [10 10];
axesFontSize = fs;
legendFontSize = fs;
% figure
% set(gcf, 'Units', 'centimeters'); 
% set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
% semilogy(Lvec(4, :), (fres_vec_poly(4, :)), '--*', 'linewidth', lw,'MarkerSize', mss)
% hold on
% semilogy(Lvec(4, :), (fres_vec_line(4, :)), '--o','linewidth', lw, 'MarkerSize', mss)
% hold on
% semilogy(Lvec(4, :), (fres_vec_FD(4, :)), '--d','linewidth', lw, 'MarkerSize', mss)
% hold on
% semilogy(Lvec(4, :), (fran_vec_poly(4, :)), '-*', 'linewidth', lw,'MarkerSize', mss)
% hold on
% semilogy(Lvec(4, :), (fran_vec_line(4, :)), '-o','linewidth', lw, 'MarkerSize', mss)
% hold on
% semilogy(Lvec(4, :), (fran_vec_FD(4, :)), '-d','linewidth', lw, 'MarkerSize', mss)
% box on
% grid on
% ylabel('Distance (m)','Interpreter', 'latex','fontsize',axesFontSize)
% xlabel('Radio stripe length (m)','Interpreter', 'latex','fontsize',axesFontSize)
% set(gca,'FontSize',fs)
% fontname(gcf, fname)
% hl2 = legend('Fresnel-Polygon','Fresnel-Line', 'Fresnel-SquareFD', ...
%     'Fraunhofer-Polygon', 'Fraunhofer-Line', 'Fraunhofer-SquareFD', 'Location', 'northwest');
% set(hl2,'interpreter','latex','FontSize',legendFontSize);

plotvec = [];
figure
subplot(1,2,1)
set(gcf, 'Units', 'centimeters'); 
set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_poly(4, :)), '-.', 'color', 'k','linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_poly(4, :)), '-', 'color', 'k', 'linewidth', lw,'MarkerSize', mss);
hold on 
plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_poly(4, :)), '-.*','color', '#0072BD', 'linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_line(4, :)), '-.o','color', '#A2142F','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_FD(4, :)), '-.d','color', '#77AC30','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_poly(4, :)), '-*', 'color', '#0072BD','linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_line(4, :)), '-o','color', '#A2142F','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_FD(4, :)), '-d','color', '#77AC30','linewidth', lw, 'MarkerSize', mss);
box on
grid on
ylabel('Distance (m)','Interpreter', 'latex','fontsize',axesFontSize)
xlabel('Radio stripe length (m)','Interpreter', 'latex','fontsize',axesFontSize)
set(gca,'FontSize',fs)
fontname(gcf, fname)
hl1 = legend([plotvec(1) plotvec(2)], '$r_{fs}$', '$r_{fr}$', 'Location', 'northwest');
ah2 = axes('position',get(gca,'position'),'visible','off');
hl2 = legend(ah2, [plotvec(6) plotvec(7) plotvec(8)], ...
    'Polygon', 'Line', 'SquareFD', 'Location', 'northwest');
set([hl1 hl2],'interpreter','latex','FontSize',legendFontSize);

% plotvec = [];
% figure
% set(gcf,'defaultAxesColorOrder',[[0.6350 0.0780 0.1840]; [0 0.4470 0.7410]]);
% set(gcf, 'Units', 'centimeters'); 
% set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
% yyaxis left
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_poly(4, :)), '*','color', 'k', 'linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_line(4, :)), 'o','color', 'k','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_FD(4, :)), 'd','color', 'k','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_poly(4, :)), '-.*', 'linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_line(4, :)), '-.o','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fran_vec_FD(4, :)), '-.d','linewidth', lw, 'MarkerSize', mss);
% hold on
% box on
% grid on
% ylabel('$r_{fr}$ (m)','Interpreter', 'latex','fontsize',axesFontSize)
% xlabel('Radio stripe length (m)','Interpreter', 'latex','fontsize',axesFontSize)
% yyaxis right
% plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_poly(4, :)), '-*','linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_line(4, :)), '-o','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Lvec(4, :), (fres_vec_FD(4, :)), '-d','linewidth', lw, 'MarkerSize', mss);
% box on
% grid on
% ylabel('$r_{fs}$ (m)','Interpreter', 'latex','fontsize',axesFontSize)
% xlabel('Radio stripe length (m)','Interpreter', 'latex','fontsize',axesFontSize)
% set(gca,'FontSize',fs)
% fontname(gcf, fname)
% hl1 = legend([plotvec(1) plotvec(2) plotvec(3)], 'Polygon', 'Line', 'SquareFD', 'Location', 'northwest');
% set(hl1,'interpreter','latex','FontSize',legendFontSize);

Fvec = 2.5e9:2.5e9:30e9;
Nvec = [9 16 25 36 49 64 81 100 121 144 169 196];

Lvec = 0.5:0.5:3;
fran_vec_poly = zeros(length(Fvec), 1);
fres_vec_poly = zeros(length(Fvec), 1);
fran_vec_line = zeros(length(Fvec), 1);
fres_vec_line = zeros(length(Fvec), 1);
fran_vec_FD = zeros(length(Fvec), 1);
fres_vec_FD = zeros(length(Fvec), 1);

Nvec = [16 36 49 64 81 100 121 144 144 169 196 196];
lambda = 3e8/2.5e9;
D1 = sqrt(2)*(sqrt(16) - 1)*(lambda/2);
for ff = 1:length(Fvec)
    f0 = Fvec(ff);
    lambda = 3e8/f0;
    inter_dist = lambda/2;
    for ll = 1:1
        maxL = 1;
        N = Nvec(ff);
        D = sqrt(2)*(sqrt(N) - 1)*inter_dist;
        fres_vec_FD(ff, 1) = ((D1^4)/(8*lambda))^(1/3);
        fran_vec_FD(ff, 1) = (2*(D1^2))/lambda;
        D = maxL;
        fres_vec_line(ff, 1) = ((D^4)/(8*lambda))^(1/3);
        fran_vec_line(ff, 1) = (2*(D^2))/lambda;
        N = floor(maxL/inter_dist);
        r0 = inter_dist/(2*sin(pi/N));
        x = zeros(N, 2);
        phii = (2*pi)/N;
        for i = 1:N
            x(i, :) = [r0*cos((i-1)*phii) r0*sin((i-1)*phii)];
        end
        dvec = zeros(N, N);
        for i = 1:N
            for j = 1:N
                dvec(i, j) = norm(x(i, :) - x(j, :));
            end
        end
        D = max(max(dvec));
        fres_vec_poly(ff, 1) = ((D^4)/(8*lambda))^(1/3);
        fran_vec_poly(ff, 1) = (2*(D^2))/lambda;
    end
end
plotvec = [];
subplot(1,2,2)
set(gcf, 'Units', 'centimeters'); 
set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_poly(:, 1)), '-.', 'color', 'k','linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_poly(:, 1)), '-', 'color', 'k', 'linewidth', lw,'MarkerSize', mss);
hold on 
plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_poly(:, 1)), '-.*','color', '#0072BD', 'linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_line(:, 1)), '-.o','color', '#A2142F','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_FD(:, 1)), '-.d','color', '#77AC30','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_poly(:, 1)), '-*', 'color', '#0072BD','linewidth', lw,'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_line(:, 1)), '-o','color', '#A2142F','linewidth', lw, 'MarkerSize', mss);
hold on
plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_FD(:, 1)), '-d','color', '#77AC30','linewidth', lw, 'MarkerSize', mss);
box on
grid on
ylabel('Distance (m)','Interpreter', 'latex','fontsize',axesFontSize)
xlabel('$f$ (GHz)','Interpreter', 'latex','fontsize',axesFontSize)
set(gca,'FontSize',fs)
fontname(gcf, fname)
hl1 = legend([plotvec(1) plotvec(2)], '$r_{fs}$', '$r_{fr}$', 'Location', 'northwest');
ah2 = axes('position',get(gca,'position'),'visible','off');
hl2 = legend(ah2, [plotvec(6) plotvec(7) plotvec(8)], ...
    'Polygon', 'Line', 'SquareFD', 'Location', 'northwest');
set([hl1 hl2],'interpreter','latex','FontSize',legendFontSize);

% plotvec = [];
% figure
% set(gcf,'defaultAxesColorOrder',[[0.6350 0.0780 0.1840]; [0 0.4470 0.7410]]);
% set(gcf, 'Units', 'centimeters'); 
% set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
% yyaxis left
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_poly(:, 1)), '*', 'color', 'k', 'linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_line(:, 1)), 'o', 'color', 'k','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_FD(:, 1)), 'd', 'color', 'k','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_poly(:, 1)), '-.*', 'linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_line(:, 1)), '-.o','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fran_vec_FD(:, 1)), '-.d','linewidth', lw, 'MarkerSize', mss);
% hold on
% box on
% grid on
% ylabel('$r_{fr}$ (m)','Interpreter', 'latex','fontsize',axesFontSize)
% xlabel('$f$ (GHz)','Interpreter', 'latex','fontsize',axesFontSize)
% yyaxis right
% plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_poly(:, 1)), '-*','linewidth', lw,'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_line(:, 1)), '-o','linewidth', lw, 'MarkerSize', mss);
% hold on
% plotvec(end + 1) = semilogy(Fvec./1e9, (fres_vec_FD(:, 1)), '-d','linewidth', lw, 'MarkerSize', mss);
% box on
% grid on
% ylabel('$r_{fs}$ (m)','Interpreter', 'latex','fontsize',axesFontSize)
% xlabel('$f$ (GHz)','Interpreter', 'latex','fontsize',axesFontSize)
% set(gca,'FontSize',fs)
% fontname(gcf, fname)
% hl1 = legend([plotvec(1) plotvec(2) plotvec(3)], 'Polygon', 'Line', 'SquareFD', 'Location', 'northwest');
% set(hl1,'interpreter','latex','FontSize',legendFontSize);