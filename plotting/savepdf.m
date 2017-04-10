function savepdf(h,fname,varargin)

figure(h);

% Use nature figure sizes: 247, 183, or 89 mm
figsize = 0;
getArgs(varargin,{'figsize=1'});

set(h,'Units','Centimeters');
pos = get(h,'Position');

set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);

set(h,'PaperPositionMode','Auto','PaperUnits','Centimeters');

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
drawPublishAxis

pos = get(h,'Position');
set(h,'Papersize',[pos(3) pos(4)]);

print(fname,'-dpdf');
