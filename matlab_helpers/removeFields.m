% removeFields.m
%
%        $Id: removeFields.m,v 1.1 2008/07/15 21:10:04 justin Exp $ 
%      usage: d = removeFields(d,removeFieldNames)
%         by: justin gardner
%       date: 07/15/08
%    purpose: removes specified fields
%             e.g.:
%             d.huh = 3;d.duh = 4;d.uhm = 5;
%             d = removeFields(d,{'huh','duh'});
%
function d = removeFields(d,removeFieldNames)

% check arguments
if ~any(nargin == [1 2])
  help removeFields
  return
end

for i = 1:length(removeFieldNames)
  if isfield(d,removeFieldNames{i})
    d = rmfield(d,removeFieldNames{i});
  end
end


