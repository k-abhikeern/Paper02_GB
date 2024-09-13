%% Need to pass on the name by which the file has to saved and it automatically saves files in a 'paper_figure' named folder
function save_fig(str)
axis square;
f = gcf;
set(f,"units","normalized");%,"outerposition",[0.0 0.0 0.8 0.8]);

cmd = ('mkdir paper_figure');
system(cmd);

str_jpg = strcat('paper_figure','/',str,'.jpg');
str_eps = strcat('paper_figure','/',str,'.eps');
str_pdf = strcat('paper_figure','/',str,'.pdf');
% exportgraphics(f,str_jpg,'Resolution',600);
% exportgraphics(f,str_eps,'Resolution',600);
% exportgraphics(f,str_pdf,'Resolution',600);

% print(f,str_jpg,'-dpng','-r600'); % Save as PNG format with 600 dpi resolution
% print(f,str_eps,'-dpng','-r600'); % Save as PNG format with 600 dpi resolution
% print(f,str_pdf,'-dpng','-r600'); % Save as PNG format with 600 dpi resolution

saveas(gcf,str_jpg);
saveas(gcf,str_pdf);
print(str_eps, '-depsc');
end