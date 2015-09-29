% calcPercentDone.m
%
%        $Id:$ 
%      usage: calcPercentDone()
%         by: justin gardner
%       date: 12/07/09
%    purpose: computes percent done for use with disppercent. For example,
%             if you have two nested loops like:
%
%             for i = 1:imax
%               for j = 1:jmax
%                 disppercent(calcPercentDone(i,imax,j,jmax);
%               end
%             end
%
function p = calcPercentDone(varargin)

% check arguments
if (nargin < 2) || ~iseven(nargin)
  help percentDone
  return
end

p = 0;
currentInterval = 1;
% loop through each set of arguments
for i = 1:2:nargin
  % add how much done we have
  p = p + currentInterval*(varargin{i}-1)/varargin{i+1};
  % the next interval will be a subdivision of this interval
  currentInterval = currentInterval/varargin{i+1};
end

  

