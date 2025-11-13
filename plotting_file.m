lw = 2.5;
fs = 14;
fname = 'Times New Roman';
mss = 8;
afFigurePosition = [5 5 13 7]; 
psize = [10 10];
axesFontSize = fs;
legendFontSize = fs;
set(gcf, 'Units', 'centimeters'); 
set(gcf, 'Position', afFigurePosition,'PaperSize',psize,'PaperPositionMode','auto');
plot(ET_vec, pow2db(Result_Fair), ':ok' ,'LineWidth', lw, 'MarkerSize', mss)
hold on 
plot(ET_vec, pow2db(Result_Cheby), '-sk', 'LineWidth', lw, 'MarkerSize', mss)
hl1 = legend('Fair-Chebyshev', 'K-Chebyshev');
box on
grid on
set(gca,'FontSize',fs)
fontsize(gcf, fs,"points")
fontname(gcf, fname)
set(hl1,'interpreter','latex','FontSize',legendFontSize);
ylabel('$\sum_{k = 1}^K 1/\Delta_{k,i}$ (dB)','Interpreter', 'latex','fontsize',axesFontSize)
xlabel('$K$','Interpreter', 'latex','fontsize',axesFontSize)
set([hl1],'color','none');
