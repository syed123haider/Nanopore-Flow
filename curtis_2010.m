clear all;
close all
clc;
%fig=figure('Position', [200 200 480 280]);

frequency = [106983 74638 24596 5168 2515 1261 500 170 25 8 3 0 0];
total_pores=max(cumsum(frequency));
percentage=frequency./total_pores;
percentage=cumsum(percentage);
percentage=flip(percentage).*100;
fig=figure('Renderer', 'painters', 'Position', [200 200 480 350]);
te=13;

left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
prices = [0 0 3 8 25 170 500 1261 2515 5168 24596 74638 106983];
yyaxis left
%bar(prices,'FaceColor',[0.9290 .6940 .1250])
bar(prices,'FaceColor',[0.5 0.5 0.5])
ylabel('$$\mathrm{Number~of~pores}$$','interpreter','latex','Color','black');
set(gca,'yscale','log','xticklabel',{'400','250','160','100','60','40','25','16','10','6','4','2.5','2'});
yyaxis right
plot(percentage,'linewidth',1.5,'color','red')
ylim([0 100]);
grid on
box on
ax = gca; % current axes
ax.FontSize = 14;
ax.TickDir = 'in';
ax.FontWeight = 'normal';
set(gca, 'XDir','reverse','ycolor','red')
xlabel('$$\mathrm{Nanopore~radius,nm}$$','interpreter','latex','Color','black')
ylabel('$$\mathrm{Percentage~of~pores}$$','interpreter','latex','Color','red')
print('-depsc2','-r400','Figures/curtis2_2010.eps');
print('-depsc2','-r400','curtis_20102.eps');

