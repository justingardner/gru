% getnum.m
%
%      usage: getnum()
%         by: justin gardner
%       date: 11/11/05
%    purpose: gets a number from the user
%             str is the text string to display to user
%             range is the valid numbers
%       e.g.: getnum('Pick a number between 1 and 3',1:3)
%

function r = getnum(str,range)

% check arguments
if ~any(nargin==[1 2])
  help getnum
  return
end

r = [];
% while we don't have a valid answer, keep asking
while(isempty(r))
  % get user input
  r = input(str,'s');
  % make sure it is a string
  if isstr(r)
    % convert to number
    r = str2num(r);
  else
    r = [];
  end
  % check if in range
  if ~isempty(r) && ~ieNotDefined('range')
    for i = 1:length(r)
      if ~any(r(i)==range)
	disp(sprintf('(getnum) %i is out of range',r(i)));
	r = [];
	break;
      end
    end
  end
end

