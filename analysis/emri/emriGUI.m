% emriGUI.m
%
%        $Id:$ 
%      usage: emriGUI()
%         by: justin gardner (based on pRFGUI)
%       date: 12/16/2022
%    purpose: GUI for getting params for emriAnal
%
function params = emriGUI(varargin)

% get the arguments
params=[];groupNum=[];defaultParams=[];scanList = [];v = [];
getArgs(varargin,{'params=[]','groupNum=[]','defaultParams=0','scanList=[]','v=[]'});

% if called with params, then just display
if ~isempty(params)
  retval = dispParams(params);
  if isempty(retval),params = [];end
  return
end

% get a view
deleteViewOnExit = false;
if isempty(v),v = newView;deleteViewOnExit = true;end

% get the group names put on top passed in group if set
groupNames = putOnTopOfList(viewGet(v,'groupName',groupNum),viewGet(v,'groupNames'));

% set the parameter string
paramsInfo = {};
paramsInfo{end+1} = {'groupName',groupNames,'Name of group from which to do pRF analysis'};
paramsInfo{end+1} = {'saveName','emriAnal','string','File name to try to save as'};

% parameters
paramsInfo{end+1} = {'frequencyAnalysis',1,'type=checkbox','Compute frequency analysis'};
paramsInfo{end+1} = {'cyclesScanMin',1,'incdec',[-1 1],'minmax',[1 inf],'Compute frequency analysis starting at this frequency','contingent=frequencyAnalysis'};
paramsInfo{end+1} = {'cyclesScanMax',5,'incdec',[-1 1],'minmax',[1 inf],'Compute frequency analysis ending at this frequency','contingent=frequencyAnalysis'};
paramsInfo{end+1} = {'detrend',{'None','Highpass','Linear','Quadratic'},'Detrend for corAnal','contingent=frequencyAnalysis'};
paramsInfo{end+1} = {'divideByMean',0,'type=checkbox','Divide by mean of time-series (not necessary for concatenation as this has already been done','contingent=frequencyAnalysis'};
paramsInfo{end+1} = {'trigonometricFunction',{'Sine','Cosine'},'Sets which function phase will be used for frequency analysis','contingent=frequencyAnalysis'};
paramsInfo{end+1} = {'temporalFiltering',0,'type=checkbox','Temporal smoothing'};
paramsInfo{end+1} = {'temporalFilter',{'Box','Gaussian'},'Type of smoothing filter','contingent=temporalFiltering'};
paramsInfo{end+1} = {'temporalFilterWidthInMs',5,'incdec',[-1 1],'minmax',[1 inf],'Size of Box smoothing filter','contingent=temporalFiltering'};
paramsInfo{end+1} = {'spatialFiltering',0,'type=checkbox','Temporal smoothing'};
paramsInfo{end+1} = {'spatialFilter',{'None','Box','Gaussian'},'Type of smoothing filter','contingent=spatialFiltering'};
paramsInfo{end+1} = {'spatialFilterWidth',0,'incdec',[-1 1],'minmax',[1 inf],'Size of Box smoothing filter','contingent=spatialFiltering'};
paramsInfo{end+1} = {'saveFilteredTSeries',0,'type=checkbox','Save the filtered time series with the same analysis overlays','contingent=frequencyAnalysis'};



% Get parameter values
if defaultParams
  params = mrParamsDefault(paramsInfo);
else
  params = mrParamsDialog(paramsInfo,'Set emriAnal parameters');
end

% if empty user hit cancel
if isempty(params)
  if deleteViewOnExit,deleteView(v);end
  return
end

% get scans
v = viewSet(v,'groupName',params.groupName);
if ~isempty(scanList)
  params.scanNum = scanList;
elseif defaultParams
  params.scanNum = 1:viewGet(v,'nScans');
else
  params.scanNum = selectScans(v);
end
if isempty(params.scanNum)
  params = [];
  if deleteViewOnExit,deleteView(v);end
  return
end

if deleteViewOnExit,deleteView(v);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just display parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = dispParams(params)

paramsInfo = {};
% grab the parameters that are indicated in paramsInfo
if isfield(params,'paramInfo')
  % get the paramsInfo
  topParamsInfo = params.paramInfo;
  % go through each one
  for i = 1:length(topParamsInfo)
    % if it exists in the params filed then add it
    if isfield(params,topParamsInfo{i}{1})
      % and it to paramInfo
      paramsInfo{end+1} = params.paramInfo{i};
      % add the value from params 
      paramsInfo{end}{2} = params.(topParamsInfo{i}{1});
      % make it non editable
      paramsInfo{end}{end+1} = 'editable=0';
    end
  end
end

retval = mrParamsDialog(paramsInfo,'pRF parameters');
retval = [];

