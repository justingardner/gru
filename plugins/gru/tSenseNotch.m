% tSenseNotch.m
%
%        $Id: tSenseNotch.m 2627 2012-10-10 07:54:28Z justin $
%      usage: tSenseNotch(v, params)
%         by: justin gardner
%       date: 11/21/2012
%    purpose: run tSense notch filter on data
%
%             to just get a default parameter structure:
% 
%             v = newView;
%             [v params] = tSenseNotch(v,[],'justGetParams=1');
%             [v params] = tSenseNotch(v,[],'justGetParams=1','defaultParams=1');
%             [v params] = tSenseNotch(v,[],'justGetParams=1','defaultParams=1','scanList=[1 2]');
%
%
function [v params] = tSenseNotch(v,params,varargin)

% check arguments
if ~any(nargin == [1 2 3 4 5 6 7 8])
  help tSenseNotch
  return
end

% other arguments
eval(evalargs(varargin));
if ieNotDefined('justGetParams'),justGetParams = 0;end
if ieNotDefined('defaultParams'),defaultParams = 0;end

% check to see if any scans have a tSense that is not one, if not 
% then warn the user that they may not need to run tSense notch
needNotchFilter = false;
defaultTSenseAcc = 1;
for iScan = 1:viewGet(v,'nScans','Raw')
  tSense = viewGet(v,'auxParam','tSense',iScan,'Raw');
  if iscell(tSense),tSense = cell2mat(tSense);end
  if isscalar(tSense) && (tSense > 1)
    needNotchFilter = true;
    defaultTSenseAcc = viewGet(v,'auxParam','tSense',iScan,'Raw');
    break;
  end
end

if ~needNotchFilter
  if ~askuser('(tSenseNotch) Did not find any scans with tSense auxParam set. Are you sure you need to run the tSense filter?')
    return
  end
end

% description of paramaters (used by mrParamsDialog functions)
paramsInfo = {...
    {'groupName',putOnTopOfList(viewGet(v,'groupName'),viewGet(v,'groupNames')),'Name of group from which scans will be taken to run tSense notch filter (usually Raw)'},...
    {'newGroupName','Notch','Name of group to which tSense Notch will be saved. If group does not already exist, it will be created.'}};
paramsInfo{end+1} = {'tSenseAcc',defaultTSenseAcc,'incdec=[-1 1]','minmax=[1 inf]','This selects what filter to use which depends on the acceleration factor has been used on the tSense data.'};
paramsInfo{end+1} = {'notchFilterHalfWidth',0,'incdec=[-1 1]','minmax=[0 inf]','Sets the halfwidth of the notch filter. If 0, then knocks out just the one frequency needed. If greater than 0, will notch out that many frequencies on either side. e.g. set to 2 will knock out 2 frequencies on either side resulting in five frequencies being knocked out altogether'};
paramsInfo{end+1} = {'dispFilter',0,'type=checkbox','Displays the filter that is being used to notch'};

    
% First get parameters
if ieNotDefined('params')
  % Initialize analysis parameters with default values
  if defaultParams
    params = mrParamsDefault(paramsInfo);
  else
    params = mrParamsDialog(paramsInfo);
  end
  % no params means user hit cancel
  if isempty(params),return,end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % select scans
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  v = viewSet(v, 'groupName', params.groupName);
  if ~ieNotDefined('scanList')
    params.scanList = scanList;
  elseif defaultParams
    params.scanList = 1:viewGet(v,'nScans');
  else
    params.scanList = selectScans(v);
  end
  if isempty(params.scanList),return,end

  % check the parameters
  params = mrParamsReconcile(params.groupName,params);
else
  % Reconcile params with current status of group and ensure that it has
  % the required fields. 
  params.paramsInfo = paramsInfo;
  params = mrParamsReconcile(params.groupName,params);
end
drawnow;

% Abort if params empty
if ieNotDefined('params'),return,end

% if just getting params then return
if justGetParams,return,end
  
% Open new view with the base group
viewBase = newView;
groupNum = viewGet(viewBase,'groupNum',params.groupName);

if (groupNum == 0)
  mrErrorDlg('tSenseNotch: ',groupName,' does not exist.');
end
viewBase = viewSet(viewBase,'currentGroup',groupNum);

% Open new view and set its group to the concat group name. Create the
% group if necessary.
viewTSeriesNotch = newView;
concatGroupNum = viewGet(viewTSeriesNotch,'groupNum',params.newGroupName);
if isempty(concatGroupNum)
  v = viewSet(v,'newgroup',params.newGroupName);
  concatGroupNum = viewGet(viewTSeriesNotch,'groupNum',params.newGroupName);
end
viewTSeriesNotch = viewSet(viewTSeriesNotch,'currentGroup',concatGroupNum);

set(viewGet(v,'figNum'),'Pointer','watch');drawnow;
tic
% Compute output volume
waitHandle = mrWaitBar(0,'(tSenseNotch) Running tSense notch filter on tSeries.  Please wait...');
v = viewSet(v,'curGroup',groupNum);
for iscan = 1:length(params.scanList)
  scanNum = params.scanList(iscan);

  viewBase = viewSet(viewBase,'curScan',params.scanList(iscan));
  v = viewSet(v,'curScan',scanNum);

  % make sure it is not a concat, which we need to write code to handle
  concatInfo = viewGet(v,'concatInfo',scanNum);
  if ~isempty(concatInfo)
    disp(sprintf('!!! (tSenseNotch) tSenseNotch is not implemented for concatenations yet. Skipping scan %s:%i !!!',viewGet(v,'groupName'),scanNum));
    continue;
  end
  
  % Load it
  mrDisp(sprintf('\n(tSenseNotch) Loading scan %i from %s\n',scanNum,viewGet(viewBase,'groupName')));
  tSeries = loadTSeries(viewBase,scanNum,'all');
	
  % run the notch filter
  [tSeries notchFilter] = tSenseNotchFilter(tSeries,params);
  
  % get the path and filename
  [path,filename,ext] = fileparts(viewGet(viewBase,'tseriesPath',scanNum));
  baseGroupName = viewGet(viewBase,'groupName');

  % Save TSeries, using scanParams from original
  scanParams = viewGet(viewBase,'scanparams');
  scanParams.fileName = [];
  scanParams.description = sprintf('Notch of %s:%i %s',baseGroupName,scanNum,viewGet(viewBase,'description',scanNum));
  scanParams.originalFileName{1} = filename;
  scanParams.originalGroupName{1} = baseGroupName;

  hdr = cbiReadNiftiHeader(viewGet(viewBase,'tseriesPath',scanNum));
  % data *MUST* be written out as float32 b/c of the small values-epm
  hdr.datatype = 16;
  
  [viewTSeriesNotch,tseriesFileName] = saveNewTSeries(viewTSeriesNotch,tSeries,scanParams,hdr);
  % get new scan number
  saveScanNum = viewGet(viewTSeriesNotch,'nScans');
    
  % Update waitbar
  mrWaitBar(iscan/length(params.scanList),waitHandle);
end
mrCloseDlg(waitHandle);
toc;

% Save evalstring for recomputing and params
evalstr = ['v = newView; v = tSenseNotch(v,params);'];
tseriesdir = viewGet(viewTSeriesNotch,'tseriesdir');
[pathstr,filename,ext] = fileparts(fullfile(tseriesdir,tseriesFileName));
save(fullfile(pathstr,filename),'evalstr','params','notchFilter');

% Delete temporary viewBase and viewTSeriesNotch
deleteView(viewBase);
deleteView(viewTSeriesNotch);

set(viewGet(v,'figNum'),'Pointer','arrow');drawnow





 
