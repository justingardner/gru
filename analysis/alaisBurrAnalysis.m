% alaisBurrAnalysis.m
%
%      usage: alaisBurrAnalysis()
%         by: justin gardner
%       date: 10/11/19
%    purpose: analyze data from Alais & Burr replication (experiment alaisburr.m)
%
function [fit e] = alaisBurrAnalysis(stimfileName,varargin)

% default return argument
fit = [];

% check arguments
if nargin < 1
  help alaisBurrAnalysis
  return
end

% parse arguments
getArgs(varargin,{'dispFit=1'});

% load and parse the stimfile
e = loadStimfile(stimfileName);
if isempty(e),return,end

% display psychometric function
fit = fitPsychometricFunction(e,dispFit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    fitPsychometricFunction    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fit = fitPsychometricFunction(e,dispFit)

% get the differences
posDiff = e.parameter.posDiff;
uniquePosDiff = unique(posDiff);

% compute whether answer are correct
correct = e.parameter.centerWhich == e.response;

% bin and average
for iVal = 1:length(uniquePosDiff)
  % get the trials with the setting of posDiff
  whichTrials = find(posDiff == uniquePosDiff(iVal));
  % compute how many trials
  nTrials = length(whichTrials);
  % compute correct
  correctBinned(iVal) = sum(correct(whichTrials))/nTrials;
  % compute ste
  correctBinnedError(iVal) = correctBinned(iVal)*(1-correctBinned(iVal))/sqrt(nTrials);
end

% fit a cumulative gaussian to data
fit = fitCumulativeGaussian(uniquePosDiff,correctBinned);

if dispFit
  % plot the figure
  mlrSmartfig('alaisBurrAnalysis');
  % plot fit
  plot(fit.fitX,fit.fitY*100,'r-');hold on
  % plot in percentile
  myerrorbar(uniquePosDiff,100*correctBinned,'yError',100*correctBinnedError,'Symbol','ko','MarkerFaceColor','k');
  xlabel('Position difference (deg)');
  ylabel('Rightward choices (percent)');
  % set title
  title(sprintf('%s (Mean: %0.2f Std: %0.2f lambda: %0.2f)',e.stimulusType,fit.mean,fit.std,fit.lambda));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%    loadStimfile    %
%%%%%%%%%%%%%%%%%%%%%%%%%
function e = loadStimfile(stimfileName)

% default return argument
e = [];

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
e = getTaskParameters(s.myscreen,s.task);
e = e{1};

% print what we found
disp(sprintf('(alaiasBurrAnalysis:loadStimfile) Found task: %s (%i trials) SID: %s Date: %s',taskFilename,e.nTrials,s.myscreen.SID,s.myscreen.starttime));

% get the variables
e.myscreen = s.myscreen;
e.task = s.task;
e.stimulus = s.stimulus;

% get experiment type
e.stimulusType = '';
if e.stimulus.visual
  e.stimulusType = sprintf('%s%s ',e.stimulusType,'Visual');
end
if e.stimulus.auditory
  e.stimulusType = sprintf('%s%s ',e.stimulusType,'Auditory');
end
if e.stimulus.bimodal
  e.stimulusType = sprintf('%s%s ',e.stimulusType,'Bimodal');
end