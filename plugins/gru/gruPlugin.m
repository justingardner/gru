% gruPlugin.m
%
%        $Id:$ 
%      usage: gruPlugin(action,<v>)
%         by: justin gardner
%       date: 12/08/2011
%    purpose: Plugin function for gru directory.
%
function retval = gruPlugin(action,v)

% check arguments
if ~any(nargin == [1 2])
  help gruPlugin
  return
end

switch action
 case {'install','i'}
  % check for a valid view
  if (nargin ~= 2) || ~isview(v)
     disp(sprintf('(gruPlugin) Need a valid view to install plugin'));
  else
    % if the view is valid, then use mlrAdjustGUI to adjust the GUI for this plugin.
    
    % this installs a new menu item called 'Select Plugins' under /Edit/ROI with the
    % separator turned on above it. It sets the callback to selectPlugins defined below.
    mlrAdjustGUI(v,'add','menu','FID Info','/Edit/Scan/Info','Callback',@callmlrFidInfo,'Separator','off');
    
    % add menu item for tSenseNotch analysis
    mlrAdjustGUI(v,'add','menu','tSense Notch','/Analysis/Motion Compensation','Callback',@callTSenseNotch,'Separator','off');
    % Install default interrogators
%    mlrAdjustGUI(v,'add','interrogator',{'gruFit'});

    % add menu item for emri analysis
    mlrAdjustGUI(v,'add','menu','emriAnal','/Analysis/Correlation Analysis','Callback',@callEmriAnal,'Separator','off');

    % add menu to analyze overlays
    mlrAdjustGUI(v,'add','menu','Analyze','/Edit/Overlay/Edit Overlay','Callback',@callOverlayAnalysis,'Separator','on');

    % This is a command that could be used to install some default colormaps
    % that will show up when you do /Edit/Overlay
    %mlrAdjustGUI(v,'add','colormap','gray');

    % This is a command that could be used to set a property of an existing menu item
    mlrAdjustGUI(v,'set','/Edit/Scan/Dicom Info','Visible','off');
    mlrAdjustGUI(v,'set','/Edit/Scan/Info','Separator','on');

    % return true to indicate successful plugin
    retval = true;
   end
 % return a help string
 case {'help','h','?'}
   retval = 'Adds UI items for handling RIKEN BSI data.';
 otherwise
   disp(sprintf('gruPlugin) Unknown command %s'));
end

%%%%%%%%%%%%%%%%%%%%%%%%
%    callmlrFidInfo    %
%%%%%%%%%%%%%%%%%%%%%%%%
function callmlrFidInfo(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% call mlrFidInfo
mlrFidInfo(v);

%%%%%%%%%%%%%%%%%%%%%%%%%
%    callTSenseNotch    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function callTSenseNotch(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% call tSenseNotch
tSenseNotch(v);

%%%%%%%%%%%%%%%%%%%%%%
%    callEmriAnal    %
%%%%%%%%%%%%%%%%%%%%%%
function callEmriAnal(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% call emriAnal
emriAnal(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    callOVerlayAnalysis    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function callOverlayAnalysis(hObject,eventdata)

% code-snippet to get the view from the hObject variable. Not needed for this callback.
v = viewGet(getfield(guidata(hObject),'viewNum'),'view');

% call overlayAnalysis
overlayAnalysis(v);

