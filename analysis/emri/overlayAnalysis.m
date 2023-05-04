% overlayAnalysis.m
%
%        $Id:$ 
%      usage: overlayAnalysis(v)
%         by: justin gardner
%       date: 01/30/23
%    purpose: Called as a plugin from mrTools Edit/Overlay menu to
%             do things like average together overlays
%
function retval = overlayAnalysis(v)

% check arguments
if ~any(nargin == [1])
  help overlayAnalysis
  return
end

% get the number of overlays
nOverlays = viewGet(v,'nOverlays');

% if there are no overlays, then nothing to do.
if nOverlays == 0
  disp(sprintf('(overlayAnalysis) No overlays loaded into viewer'));
  return
end

% put up dialog box with parameter choices. 
paramsInfo = [];
paramsInfo{end+1} = {'analysisType',{'Average'},'Average selected overlays'};
params = mrParamsDialog(paramsInfo,'Select analysis parameters');
if isempty(params)
  disp(sprintf('(overlayAnalysis) Analysis aborted'));
  return
end

% now choose overlays
overlayList = selectInList(v,'overlay',sprintf('Select overlays for: %s',params.analysisType));
if isempty(overlayList)
  disp(sprintf('(overlayAnalysis) No overlays selected. Aborting.'));
  return
end

% load the overlays
inputOverlays = viewGet(v,'overlay',overlayList);
if isempty(inputOverlays)
  disp(sprintf('(overlayAnalysis) Overlays %s not found. Aborting.',mlrnum2str(overlayList)));
  return
end

% run the analysi
switch (params.analysisType)
  case {'Average'}
   outputOverlay = averageOverlays(inputOverlays);
 otherwise
  disp(sprintf('(overlayAnalysis) Analysis %s not known',params.analysisType));
end

% checkn for empty analysis
if isempty(outputOverlay)
  disp(sprintf('(overlayAnalysis) Empty overlay returned for %s. Aborting',params.analysisType))
  return
end

% install the overlays
for iOverlay = 1:length(outputOverlay)
  % set standard fields for the overlays
  outputOverlay(iOverlay).function = params.analysisType;
  v = viewSet(v,'newOverlay',outputOverlay(iOverlay));
end
  


%%%%%%%%%%%%%%%%%%%%%%%%%
%    averageOverlays    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function outputOverlay = averageOverlays(inputOverlays)

outputOverlay = [];

% check to make sure there are enough overlays
if length(inputOverlays) < 2
  disp(sprintf('(overlayAnalysis:averageOverlays) Need more than %i overlays to average',length(inputOverlays)));
  return
end

% initialize outputOverlay with one of the input overlays
outputOverlay = inputOverlay(1);

% TODO: average overlays into the outputOverlay and return
