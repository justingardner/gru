% gruSubjectnum2ID
%
%      usage: id = gruSubjectNum2ID(num)
%         by: justin gardner
%       date: 04/27/15
%    purpose: convert subject num like 25 to id: s0025
%
function id = gruSubjectNum2ID(num)

id = [];

% check arguments
if ~any(nargin == [1])
  help gruSubjectID2n
  return
end

% already is a string
if isstr(num)
  id = num;
  return
end

% convert it
id = sprintf('s%04i',num);


