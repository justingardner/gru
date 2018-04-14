function savepdf(h,fname)

figure(h);

% Use nature figure sizes: 247, 183, or 89 mm

% For 4-figure 247 use 5.8 width and 10 height
% For 3-figure 247 use 8 width and 10 height

l = legend(gca,'boxoff');
set(l,'Color','none');
            
set(h,'Units','Centimeters');
pos = get(h,'Position');

set(gcf,'Color',[1 1 1]);
set(gca,'Color',[1 1 1]);

set(h,'PaperPositionMode','Auto','PaperUnits','Centimeters');

set(h,'Papersize',[pos(3) pos(4)]);

print(fname,'-dpdf');
