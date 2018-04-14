% dline.m
%
%      usage: dline(<c>,<oneToOne=0>)
%         by: justin gardner
%       date: 08/22/03
%       e.g.: dline;
%    purpose: plot a diagonal line, c is an optional color/symbol e.g. 'k-'
%             if oneToOne is set to 1, then the 1-1 diagonal will be plotted
%             otherwise the diagonal for the current axis settings will be plotted (i.e.
%             from the bottom left to top right corner
%
function retval = dline(c,oneToOne)

% check arguments
if ~any(nargin == [0 1 2])
  help dline
  return
end

if ieNotDefined('c'),c = 'k:';end
if ieNotDefined('oneToOne');oneToOne = 0;end

ax = axis(gca);
minx = ax(1);maxx = ax(2);
if isequal(get(gca,'XScale'),'log')
  minx = min(get(gca,'Xtick'));
end

miny = ax(3);maxy = ax(4);
if isequal(get(gca,'YScale'),'log')
  miny = min(get(gca,'Ytick'));
end

if oneToOne
  retval = plot([max(minx,miny) min(maxx,maxy)],[max(minx,miny) min(maxx,maxy)],c);
else
  retval = plot([minx maxx],[miny maxy],c);
end  