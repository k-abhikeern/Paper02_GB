%% Need to pass legend name in the main file(if not passed anything then leave it blank i.e.setFigureProperties())
function setFigureProperties2(lgd)
    set(groot, 'DefaultAxesFontSize', 16); % Set default font size
    set(groot, 'DefaultAxesFontName', 'SansSerif'); % Set default font style
    set(groot, 'DefaultTextFontSize', 16); % Set default font size for text elements
    set(groot, 'DefaultTextFontName', 'SansSerif'); % Set default font style for text elements
    set(groot, 'DefaultLegendFontSize', 16); % Set default font size for legend
    set(groot, 'DefaultLegendFontName', 'SansSerif'); % Set default font style for legend
    set(groot, 'DefaultTextInterpreter', 'latex'); % Set interpreter to LaTeX
    set(gca,'fontsize', 16);
    
%% Legend settings
% Copy following lines to the the main file just below where you want to set legend
%     lgd = legend('NEMD Data');
%     setFigureProperties(lgd);

    lgd.FontName = 'SansSerif';
    lgd.Interpreter = 'latex';
    lgd.FontSize = 16;
    lgd.Location = 'northwest';
    lgd.NumColumns = 1;
end
