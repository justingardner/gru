% splitcmap.m
%
%      usage: cmap = splitcmap(<n>)
%         by: justin gardner
%       date: 04/20/07
%    purpose: return a split colormap
%
function cmap = splitcmap(n)

% check arguments
if ~any(nargin == [0 1])
  help splitcmap
  return
end

% make sure n is defined properly
if ieNotDefined('n'),n=256;end
n = round(n);

% create a split colorbar
% take the upper part of the hot map
% the number 140 here is arbitrarily chosen
% to make the map start at a good color
hotcmap = hot(round((140/128)*(n/2)));
hotcmap = hotcmap(end-(floor(n/2)-1):end,:);

% flip the map for the other side
coolcmap(:,1) = hotcmap(:,3);
coolcmap(:,2) = hotcmap(:,2);
coolcmap(:,3) = hotcmap(:,1);

% now put the two together, adding an extra
% element in the middle of an odd number
% many colors, so fix that here
if isodd(n)
  cmap = [flipud(coolcmap);[0 0 0];hotcmap];
else
  cmap = [flipud(coolcmap);hotcmap];
end
  

