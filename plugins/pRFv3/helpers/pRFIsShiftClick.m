% pRFIsShiftClick
%
%      usage: shiftDown = pRFIsShiftClick(selectionType,currentModifier)
%         by: Austin Kuo with Codex
%       date: 09/02/2026
%    purpose: Combine the two modifier signals available on MATLAB figures.
%
function shiftDown = pRFIsShiftClick(selectionType,currentModifier)

if nargin < 1 || isempty(selectionType)
  selectionType = '';
end
if nargin < 2 || isempty(currentModifier)
  currentModifier = {};
end

modifierMatch = strcmpi(currentModifier,'shift');
shiftDown = strcmpi(selectionType,'extend') || any(modifierMatch(:));
