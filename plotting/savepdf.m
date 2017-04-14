function savepdf(h,fname,varargin)

figure(h);

% Use nature figure sizes: 247, 183, or 89 mm

% For 4-figure 247 use 5.8 width and 10 height
% For 3-figure 247 use 8 width and 10 height


%% Change the figure size
figsize = 0; subplots = 0;
getArgs(varargin,{'figsize=0','subplots=0'});

set(h,'Units','Centimeters');
pos = get(h,'Position');

set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);

set(h,'PaperPositionMode','Auto','PaperUnits','Centimeters');

if length(figsize)==2
    pos(4) = figsize(2);
    pos(3) = figsize(1);
    figsize = figsize(1);
end

switch figsize
    case 0
        % do nothing
        figsize = pos(3);
    case 0.5
        % half column
        figsize = 8.9;
    case 1
        % full column
        figsize = 18.3;
    case 2 
        % double column
        figsize = 24.7;
end

set(h,'Position',[0, 0, figsize, pos(4)*figsize/pos(3)]);

% DPA
if ~subplots
    % Remove the legend box + background
    l = legend(gca,'boxoff');
    set(l,'Color','none');
    drawPublishAxis('calledbysavepdf=1');
else
    for x = 1:subplots(1)
        for y = 1:subplots(2)
            subplot(subplots(1),subplots(2),(x-1)*subplots(2)+y);
            % Remove the legend box + background
            l = legend(gca,'boxoff');
            set(l,'Color','none');
            drawPublishAxis('calledbysavepdf=1');
        end
    end
end

pos = get(h,'Position');
set(h,'Papersize',[pos(3) pos(4)]);

print(fname,'-dpdf');
