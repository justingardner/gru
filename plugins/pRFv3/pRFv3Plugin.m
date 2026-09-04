% pRFv3Plugin.m
%
%        $Id:$
%      usage: pRFv3Plugin(action,<v>)
%         by: justin gardner
%       date: 11/24/10
%    purpose: Plugin function for pRF directory.
%
function retval = pRFv3Plugin(action,v)

retval = false;

% check arguments
if ~any(nargin == [1 2])
  help pRFv3Plugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(pRFv3Plugin) Need a valid view to install plugin'));
  else
    % pRFv2 and pRFv3 expose the same public MATLAB function names and must
    % never be installed together. If both were checked, prefer the newer
    % plugin, persist that correction, and hide the v2 menu that mlrPlugin
    % may already have installed earlier in this view.
    if makePRFv3SelectionExclusive
      mlrAdjustGUI(v,'remove','menu','pRF Analysis (v2)');
    end

    % pRFv2 and pRFv3 intentionally expose the same analysis function
    % names. Selecting pRFv3 makes its runtime paths take precedence.
    activatePRFv3Path;

    % Replace listeners from an older pRF version or a prior plugin load,
    % then guard the legacy R2025b getpts callback on this mrLoadRet figure.
    pRFInstallGetptsCompatibility('teardown');
    pRFInstallGetptsCompatibility('install',viewGet(v,'figureNumber'));

    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.

    % this installs a new menu item called 'Select Plugins' under /Edit/ROI with the
    % separator turned on above it. It sets the callback to selectPlugins defined below.
    mlrAdjustGUI(v,'add','menu','pRF Analysis (v3)','/Analysis/Correlation Analysis','Callback',@callpRF,'Separator','off');

    % Install default interrogators
    mlrAdjustGUI(v,'add','interrogator',{'pRFFit'});

    % This is a command that could be used to install some default colormaps
    % that will show up when you do /Edit/Overlay
    %mlrAdjustGUI(v,'add','colormap','gray');

    % This is a command that could be used to set a property of an existing menu item
    %mlrAdjustGUI(v,'set','Plots/Mean Time Series','Separator','on');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Runs population receptive field analysis (mutually exclusive with pRFv2).';
 case {'uninstall','u','remove'}
   % mlrPlugin does not currently issue an uninstall callback, but expose a
   % cleanup action for explicit plugin reloads and path switching.
   activatePRFv3Path;
   pRFInstallGetptsCompatibility('teardown');
   retval = true;
 otherwise
   disp(sprintf('(pRFv3Plugin) Unknown command %s',action));
end

%%%%%%%%%%%%%%%%%%%%%%%
%    callpRF    %
%%%%%%%%%%%%%%%%%%%%%%%
function callpRF(hObject,eventdata)

activatePRFv3Path;

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');
v = pRF(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    activatePRFv3Path    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pluginRoot = activatePRFv3Path

% mlrPlugin only requires the plugin entry point itself to be on the path.
% Bootstrap the tree relative to that entry point before calling a helper
% that may live in any pRFv3 subdirectory.
pluginRoot = fileparts(mfilename('fullpath'));
pluginTree = genpath(pluginRoot);
if isempty(pluginTree)
  error('pRFv3Plugin:MissingPluginTree', ...
    'Could not construct the MATLAB path for %s.',pluginRoot);
end
addpath(pluginTree,'-begin');
pluginRoot = pRFv3ActivatePath(pluginRoot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    makePRFv3SelectionExclusive    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changedSelection = makePRFv3SelectionExclusive

changedSelection = false;
try
  selectedPlugins = mrGetPref('selectedPlugins');
catch preferenceException
  error('pRFv3Plugin:CannotReadSelectedPlugins', ...
    'Could not verify pRF plugin exclusivity: %s',preferenceException.message);
end

if ischar(selectedPlugins)
  if isempty(selectedPlugins)
    selectedPlugins = {};
  else
    selectedPlugins = regexp(selectedPlugins,',','split');
  end
elseif isstring(selectedPlugins)
  selectedPlugins = cellstr(selectedPlugins(:)');
elseif ~iscell(selectedPlugins)
  error('pRFv3Plugin:InvalidSelectedPlugins', ...
    'Could not verify pRF plugin exclusivity because selectedPlugins has an unexpected type.');
end

hasPRFv2 = any(strcmp('pRFv2',selectedPlugins));
if ~hasPRFv2
  return
end

selectedPlugins = selectedPlugins(~strcmp('pRFv2',selectedPlugins));
if ~any(strcmp('pRFv3',selectedPlugins))
  selectedPlugins{end+1} = 'pRFv3';
end
try
  mrSetPref('selectedPlugins',selectedPlugins,false);
catch preferenceException
  error('pRFv3Plugin:CannotEnforceMutualExclusion', ...
    ['pRFv2 and pRFv3 cannot be loaded together, and pRFv2 could not be ' ...
     'removed from selectedPlugins: %s'],preferenceException.message);
end

changedSelection = true;
disp('(pRFv3Plugin) pRFv2 and pRFv3 are mutually exclusive.');
disp('(pRFv3Plugin) Removed pRFv2 from selectedPlugins; pRFv3 will be used.');
