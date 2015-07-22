% gruSubjectID2num
%
%        $Id:$ 
%      usage: num = gruSubjectID2n(id)
%         by: justin gardner
%       date: 04/27/15
%    purpose: conver subject id like s0001 into number (ie 1)
%
function num = gruSubjectID2num(id)

num = [];

% check arguments
if ~any(nargin == [1])
  help gruSubjectID2n
  return
end

% already is a number
if ~isstr(id)
  num = id;
  return
end

% check the first letter is s
if (length(id) > 1) && isequal(lower(id(1)),'s')
  num = str2num(id(2:end));
end

