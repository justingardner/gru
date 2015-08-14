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
%             dialogs to choose scans that you want to run on and
%             set parameters, then should just run.
%
function retval = mlrRunMultiPRF()

% check arguments
if ~any(nargin == [0])
  help mlrRunMultiPRF
  return
end

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
  disp(sprintf('(mlrRunMultPRF) Load the ROIs that you want to run the pRF analysis'));
  v = loadROI(v);

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
  sessionInfo(end).params = params;
end

% show user what we are about to run
dispHeader
for iSession = 1:length(sessionInfo)
  disp(sprintf('%i: %s',iSession,sessionInfo(iSession).pRFSessionPath));
end
dispHeader

% then run them
if askuser('Run the pRF Sessions listed in the command window?',0,1)
  for iSession = 1:length(sessionInfo)
    cd(sessionInfo(iSession).pRFSessionPath);
    v = newView;
    v = viewSet(v,'curGroup',sessionInfo(iSession).groupName);
    v = loadROI(v,sessionInfo(iSession).roiNames);
    v = pRF(v,sessionInfo(iSession).params);
    mrQuit;
  end
end
