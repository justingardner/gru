function activated = pRFv3Startup
% pRFv3Startup
%
%      usage: activated = pRFv3Startup
%         by: Austin Kuo with Codex
%    purpose: Honor the saved mlrPlugin selection during MATLAB startup,
%             including sessions in which no mrLoadRet window is opened.

activated = false;

if exist('mrGetPref','file') ~= 2
  error('pRFv3Startup:MissingMrGetPref', ...
    'Could not inspect selectedPlugins because mrGetPref is not on the MATLAB path.');
end

try
  selectedPlugins = mrGetPref('selectedPlugins');
catch preferenceException
  error('pRFv3Startup:CannotReadSelectedPlugins', ...
    'Could not inspect selectedPlugins: %s',preferenceException.message);
end

if isempty(selectedPlugins)
  return
elseif ischar(selectedPlugins)
  selectedPlugins = regexp(selectedPlugins,',','split');
elseif isstring(selectedPlugins)
  selectedPlugins = cellstr(selectedPlugins(:)');
elseif ~iscell(selectedPlugins)
  error('pRFv3Startup:InvalidSelectedPlugins', ...
    'selectedPlugins must be a character vector, string array, or cell array of names.');
end

validPluginNames = cellfun(@(pluginName) ischar(pluginName) || ...
  (isstring(pluginName) && isscalar(pluginName)),selectedPlugins);
if ~all(validPluginNames)
  error('pRFv3Startup:InvalidSelectedPlugins', ...
    'Every selectedPlugins entry must be a character vector or scalar string.');
end
selectedPlugins = cellfun(@char,selectedPlugins,'UniformOutput',false);
selectedPlugins = strtrim(selectedPlugins);

if ~any(strcmp('pRFv3',selectedPlugins))
  return
end

if any(strcmp('pRFv2',selectedPlugins))
  warning('pRFv3Startup:ConflictingSelection', ...
    'Both pRFv2 and pRFv3 are selected; activating pRFv3.');
end

% Let the activation helper locate the plugin root from its own position.
% This remains correct if either helper moves within the pRFv3 tree.
pRFv3ActivatePath;
activated = true;
