% mlrRunMultiPRF.m
%
%        $Id:$ 
%      usage: mlrRunMultiPRF()
%         by: justin gardner
%       date: 08/14/15
%    purpose: Run multiple pRF runs one after the other. This is
%             useful if you are processing multiple subjects data
%             but don't want to run all of them concurrently (which
%             might eat up data). This will step you through a few
%             dialogs to choose scans that you want to run on, load
%             bases to restrict on and set parameters, then should just run.
% 
%             If you want to restrict on ROIs instead of bases
%             mlrRunMultiPRF('useROIs=1');
%
%
function retval = mlrRunMultiPRF(varargin)

% arugments
getArgs(varargin,{'useROIs=0'});

% make sure we are running mrTools
mlrPath('mrTools');

% get current path
curpwd = pwd;

% set directory to start with where data lives
dataDir = fullfile(mlrReplaceTilde('~/data'),'mlrAnatDB');
if ~isdir(dataDir)
  dataDir = mlrReplaceTilde('~/data');
end

% variable that will contain session info
sessionInfo = [];

mrQuit;

% loop until the user does not have any sessions to run anymore
getSessions = true;firstLoop = true;
while (getSessions)
  if ~firstLoop
    % ask the user if this is it
    getSessions = askuser('Add another session to do pRF on?',0,1);
    if ~getSessions,continue,end
  end
  firstLoop = false;

  % get the session
  pRFSession = getPathStrDialog(dataDir,'Find pRF session to analyze (click on mrSession.mat file)','*.mat');
  if isempty(pRFSession),continue,end

  % get the path for the next round
  dataDir = fileparts(pRFSession);
  
  % cd to that session
  pRFSessionPath = fileparts(pRFSession);
  if isdir(pRFSessionPath)
    cd(pRFSessionPath);
  else
    continue
  end

  % open up a view to that session
  v = newView;
  if isempty(v),mrQuit;continue,end

  % load the ROI
  if useROIs
    disp(sprintf('(mlrRunMultPRF) Load the ROIs that you want to run the pRF analysis on (ok to load none)'));
    v = loadROI(v);
  else
    % load any anatomies
    disp(sprintf('(mlrRunMultPRF) Load any base anatomies that you want to run the pRF analysis on (ok to load none)'));
    v = loadAnat(v);
  end
  
  % select group
  groupNames = viewGet(v,'groupNames');
  paramsInfo = {{'groupName',fliplr(groupNames),'Select which group you are going to run pRF on'}};
  params = mrParamsDialog(paramsInfo,'Select which group you are going to run pRF on');
  if isempty(params),mrQuit;continue,end
  
  % change the group
  v = viewSet(v,'curGroup',params.groupName);
  
  % set the parameters
  disp(sprintf('(mlrRunMultPRF) Set parameters of pRF analysis'));
  [v params] = pRF(v,[],'justGetParams=1');
  if isempty(params),mrQuit;continue,end
  
  % keep the info
  sessionInfo(end+1).pRFSessionPath = pRFSessionPath;
  sessionInfo(end).groupName = params.groupName;
  sessionInfo(end).roiNames = viewGet(v,'roiNames');
  sessionInfo(end).baseNames = viewGet(v,'baseNames');
  sessionInfo(end).params = params;
  
  % quit the session
  mrQuit;
end

% show user what we are about to run
dispHeader
for iSession = 1:length(sessionInfo)
  % get base name string
  if strcmp(sessionInfo(iSession).params.restrict,'Base: ALL')
    % get all base names
    baseNames = '';
    for iBase = 1:length(sessionInfo(iSession).baseNames)
      baseNames = sprintf('%s%s, ',baseNames,sessionInfo(iSession).baseNames{iBase});
    end
    baseNames = baseNames(1:end-2);
    % display what we are doing with list of base names
    sessionInfo(iSession).dispStr = sprintf('%i) %s group:%s scans:%s restrict:%s (%s)',iSession,sessionInfo(iSession).pRFSessionPath,sessionInfo(iSession).groupName,num2str(sessionInfo(iSession).params.scanNum),sessionInfo(iSession).params.restrict,baseNames);
  else
    % set up display string
    sessionInfo(iSession).dispStr = sprintf('%i) %s group:%s scans:%s restrict:%s',iSession,sessionInfo(iSession).pRFSessionPath,sessionInfo(iSession).groupName,num2str(sessionInfo(iSession).params.scanNum),sessionInfo(iSession).params.restrict);
  end
  % display what we are doing
  disp(sessionInfo(iSession).dispStr);
end
dispHeader

% then run them
if askuser('Run the pRF Sessions listed in the command window?',0,1)
  for iSession = 1:length(sessionInfo)
    dispHeader
    disp(sessionInfo(iSession).dispStr);
    dispHeader
    cd(sessionInfo(iSession).pRFSessionPath);
    v = newView;
    v = viewSet(v,'curGroup',sessionInfo(iSession).groupName);
    % load any rois
    if ~isempty(sessionInfo(iSession).roiNames)
      v = loadROI(v,sessionInfo(iSession).roiNames);
    end
    % load any bases
    if ~isempty(sessionInfo(iSession).baseNames)
      v = loadAnat(v,sessionInfo(iSession).baseNames);
    end
    v = pRF(v,sessionInfo(iSession).params);
    mrQuit;
  end
end
