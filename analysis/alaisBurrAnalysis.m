% alaisBurrAnalysis.m
%
%      usage: alaisBurrAnalysis()
%         by: justin gardner
%       date: 10/11/19
%    purpose: analyze data from Alais & Burr replication (experiment alaisburr.m)
%
function e = alaisBurrAnalysis(stimfileNames,varargin)

% default return argument
fit = [];

% default to working on current directory
if nargin < 1, stimfileNames = [];end

% parse arguments
getArgs(varargin,{'dispFit=1'});

% get filenames and path
[e.path stimfileNames] = getStimfileNames(stimfileNames);
if isempty(e.path),return,end

% check for valid files as we go through
% nFiles will contain how many vaild files we have
e.nFiles = 0;

% cycle through all files
for iFile = 1:length(stimfileNames)

  % load and parse the stimfile
  d = loadStimfile(fullfile(e.path,stimfileNames{iFile}));

  % valid file, so keep its information
  if ~isempty(d)
    % update count
    e.nFiles = e.nFiles + 1;
    % and list of filenames
    e.filenames{e.nFiles} = stimfileNames{iFile};
    % and the data
    e.d{e.nFiles} = d;

    % fit psychometric function
    e.d{e.nFiles} = fitPsychometricFunction(e.d{end});
  end
end

% if no valid files found return
if e.nFiles == 0,return,end

% convert d into an array to make it easier to work with
e.d = cell2mat(e.d);

% display the fits
if dispFit
  % open figure
  mlrSmartfig('alaisBurrAnalysis_psychometricfits','reuse');clf;
  % plot each function
  for iFile = 1:e.nFiles
    % set subplot
    subplot(1,e.nFiles,iFile);
    %diplay fit
    dispPsychometricFunction(e.d(iFile));
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display psychometric function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dispPsychometricFunction(d)

% set colors to display in
if (d.nCond == 1)
  dataColors = {'k'};
  fitColors = {'r'};
else
  for iCond = 1:d.nCond
    dataColors{iCond} = getSmoothColor(iCond,d.nCond,'cool');
    fitColors{iCond} = getSmoothColor(iCond,d.nCond,'cool');
  end
end

% title string
titleStr = d.experimentName;

for iCond = 1:d.nCond
  % plot fit
  plot(d.fit(iCond).fitX,d.fit(iCond).fitY*100,'-','Color',fitColors{iCond});hold on
  
  % plot in percentile
  myerrorbar(d.cond(iCond).uniquePosDiff,100*d.cond(iCond).correctBinned,'yError',100*d.cond(iCond).correctBinnedError,'Symbol','o','MarkerFaceColor',dataColors{iCond});
  xlabel('Position difference (deg)');

  % append fit parameters to title
  titleStr = sprintf('%s\nMean: %0.2f Std: %0.2f lambda: %0.2f',titleStr,d.fit(iCond).mean,d.fit(iCond).std,d.fit(iCond).lambda);
end

% display title
title(titleStr);

% display a legend for more than one
if d.nCond > 1
  % display legend
  hLegend = mylegend(d.condNames,dataColors);
  set(hLegend,'Location','northwest');
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitPsychometricFunction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = fitPsychometricFunction(d)

% fit for each set of data in the d structure
for iCond = 1:d.nCond
  % get these trial nums
  trialNums = d.condTrialNums{iCond};

  % get the differences
  d.cond(iCond).posDiff = d.parameter.posDiff(trialNums);
  d.cond(iCond).uniquePosDiff = unique(d.cond(iCond).posDiff);

  % compute whether answer are correct
  correct = d.parameter.centerWhich(trialNums) == d.response(trialNums);

  % bin and average
  for iVal = 1:length(d.cond(iCond).uniquePosDiff)
    % get the trials with the setting of posDiff
    whichTrials = find(d.cond(iCond).posDiff == d.cond(iCond).uniquePosDiff(iVal));
    % compute how many trials
    nTrials = length(whichTrials);
    % compute correct
    d.cond(iCond).correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
    % compute ste
    d.cond(iCond).correctBinnedError(iVal) = d.cond(iCond).correctBinned(iVal)*(1-d.cond(iCond).correctBinned(iVal))/sqrt(nTrials);
  end

  % fit a cumulative gaussian to data
  d.fit(iCond) = fitCumulativeGaussian(d.cond(iCond).uniquePosDiff,d.cond(iCond).correctBinned);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function d = loadStimfile(stimfileName)

% default to empty
d = [];

% add mat extension
stimfileName = setext(stimfileName,'mat');

% see if it exists
if ~isfile(stimfileName)
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) Could not find stimfile: %s',stimfileName));
  return
end

% load the file
s = load(stimfileName);
if ~isfield(s,'myscreen') || ~isfield(s,'task')
  disp(sprintf('(alaisBurrAnalysis:loadStimfile) No myscreen or task in stimfile: %s',stimfileName));
  return
end

% now check to see if this has the alaisBurr experiment in it
if ~iscell(s.task) || ~iscell(s.task{1})
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Task variable has incorrect phases in stimfile: %s',stimfileName));
  return
end

% check task filename
taskFilename = s.task{1}{1}.taskFilename;
if isempty(strfind(lower(taskFilename),'alaisburr'))
  disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Incorrect task in stimfile: %s',taskFilename));
  return
end

% parse into parameters
d = getTaskParameters(s.myscreen,s.task);
d = d{1};

% print what we found
disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,d.nTrials,s.myscreen.SID,s.myscreen.starttime));

% get the variables
d.myscreen = s.myscreen;
d.task = s.task;
d.stimulus = s.stimulus;

% get experiment type
d.stimulusType = '';
if d.stimulus.visual
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Visual');
end
if d.stimulus.auditory
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Auditory');
end
if d.stimulus.bimodal
  d.stimulusType = sprintf('%s%s ',d.stimulusType,'Bimodal');
end
d.stimulusType = strtrim(d.stimulusType);

% get width
d.visualWidth = d.stimulus.width;
if isfield(s.task{1}{1}.parameter,'displacement')
  d.displacement = s.task{1}{1}.parameter.displacement;
else
  d.displacement = [];
end

% set experiment name and conditions
if d.stimulus.bimodal
  % set experiment name
  d.experimentName = sprintf('%s: [%s] width: %0.1f',d.stimulusType,mlrnum2str(d.displacement),d.visualWidth);
  % set number of conditions
  d.nCond = length(d.displacement);
  % get trials for each condition
  for iCond = 1:d.nCond
    % set condition names
    d.condNames{iCond} = sprintf('Displacement: %s',mlrnum2str(d.displacement(iCond)));
    % get trials for this condition
    d.condTrialNums{iCond} = find(d.parameter.displacement == d.displacement(iCond));
  end
elseif d.stimulus.visual
  % set experiment name
  d.experimentName = sprintf('%s: width: %0.1f',d.stimulusType,d.visualWidth);
  % set number of conditions
  d.nCond = 1;
  d.condNames{1} = sprintf('Visual: width: %0.1f',d.visualWidth);
  % set trial nums to all
  d.condTrialNums{1} = 1:d.nTrials;
else
  % set experiment name
  d.experimentName = sprintf('%s',d.stimulusType);
  % set number of conditions
  d.nCond = 1;
  % set condition name
  d.condNames{1} = sprintf('Auditory');
  % set trial nums to all
  d.condTrialNums{1} = 1:d.nTrials;
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%
%    getStimfileNames    %
%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stimfilePath stimfileNames] = getStimfileNames(stimfileNames)


%  if passed in empty then it means to run on current directory
if isempty(stimfileNames)
  % get current directory
  stimfilePath = pwd;
  % get everything in the directory that is a mat file
  matfiles = dir('*.mat');
  stimfileNames = {matfiles(:).name};
elseif isfile(setext(stimfileNames,'mat'))
  % make sure the extension is .mat
  stimfileNames = setext(stimfileNames,'mat');
  % first check to see if it has a path
  if ~isempty(fileparts(stimfileNames))
    stimfilePath = fileparts(stimfileNames);
    stimfileNames = getLastDir(stimfileNames);
  else
    stimfilePath = pwd;
  end
else
  disp(sprintf('(alaisBurrAnalysis) Could not find %s',stimfileNames));
  stimfilePath = '';
  stimfileNames = '';
end

% make sure we are returning a cell array
stimfileNames = cellArray(stimfileNames);


  
  
  

